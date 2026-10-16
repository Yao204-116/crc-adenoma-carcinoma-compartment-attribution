# -*- coding: utf-8 -*-
"""一条命令走完：GitHub 建仓/推送 → 打 tag → 发 Release → 取 Zenodo DOI → 回填稿件。

为什么单独写一个脚本：
  * `git push` 需要凭据。凭据**不写入** `.git/config`（用一次性带 token 的 URL 推送），
    且所有回显都做 token 脱敏 —— git 报错时会把 URL 原样打出来，不脱敏就会泄漏 PAT。
  * Zenodo 的 GitHub 集成只认 release，所以 DOI 只能在 push + release 之后轮询取回。
  * 取回 DOI 后写进 `code_availability.json`（Code availability 的单一数据源）并重建
    投稿稿，稿件里的三处占位就不再是占位。

用法（要在能出网的环境下运行）：
    python push_github.py --token-file _push/.ghtoken

    # github.com:443 不通、但 api.github.com 正常时（国内网络常见）：
    python push_github.py --token-file _push/.ghtoken --via-api

常用选项：
    --repo-name NAME        仓库名（默认 crc-adenoma-carcinoma-compartment-attribution）
    --r-version "R 4.5.1"   填进 Code availability 的 R 版本
    --license MIT
    --via-api               用 REST API 推内容与 tag（绕开被阻断的 github.com）
    --no-create             仓库已手动建好，不要尝试自动创建
    --allow-placeholder     跳过发布前闸门 preflight_release.py（默认拒绝，强烈不建议）
    --skip-release          只 push，不发 release、不取 DOI
    --poll 20               等待 Zenodo 归档的轮数（每轮 30 秒）；0 = 不等待
"""
import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.github.com"
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.join(HERE, "repo")
DEFAULT_REPO = "crc-adenoma-carcinoma-compartment-attribution"
ZENODO_API = "https://zenodo.org/api/records"

_TOKEN = ""            # 全局，供脱敏使用


def redact(s):
    """把 token 从任何要打印的文本里抹掉。"""
    if not s:
        return s
    out = str(s)
    if _TOKEN:
        out = out.replace(_TOKEN, "***REDACTED***")
        # 也防 urlencode / 大小写变体
        out = out.replace(urllib.parse.quote(_TOKEN, safe=""), "***REDACTED***")
    return out


def log(step, msg):
    print(f"[{step}] {msg}", flush=True)


def gh(method, path, payload=None, token=None):
    """调 GitHub REST；返回 (status, dict)。HTTP 错误也返回而不抛。"""
    url = path if path.startswith("http") else API + path
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {token or _TOKEN}")
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "reproducible-code-deposit")
    if data:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            body = r.read().decode("utf-8", "replace")
            return r.status, (json.loads(body) if body.strip() else {})
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            return e.code, json.loads(body)
        except Exception:
            return e.code, {"message": body[:400]}
    except Exception as e:                      # 网络层面的失败
        return 0, {"message": f"{type(e).__name__}: {e}"}


def git(*args, check=True):
    p = subprocess.run(["git", *args], cwd=REPO_DIR,
                       capture_output=True, text=True)
    if check and p.returncode != 0:
        raise SystemExit("\n".join([
            f"git {' '.join(a for a in args if _TOKEN not in a)} 失败",
            redact(p.stdout), redact(p.stderr)]))
    return p


# 会**公开发布**的文本资产。占位符必须同时填进这两个文件：
# ★ 只填 README 会漏掉 RELEASE_NOTES_v1.0.0.md —— 而 release notes 的正文正是
#   贴进 GitHub release 页面、随即被 Zenodo 归档的那段文字。实测漏网过一次，
#   是 preflight_release.py 扫全仓库占位符时才发现的。
PUBLISHED_TEXT = ["README.md", "RELEASE_NOTES_v1.0.0.md"]


def fill_placeholders(sub, commit_msg, step):
    """把 {KEY: value} 填进 PUBLISHED_TEXT 里所有 `<KEY>` / `\\<KEY\\>` 两种写法。

    markdown 里占位符写作反斜杠转义形式，故两种都要替换（只替换一种会静默漏掉
    表格行 —— 实测 README 的 "Source repository" 一格就这样漏过）。
    """
    touched = []
    for rel in PUBLISHED_TEXT:
        fp = os.path.join(REPO_DIR, rel)
        if not os.path.exists(fp):
            continue
        txt = open(fp, encoding="utf-8").read()
        new = txt
        for k, v in sub.items():
            new = new.replace("\\<%s\\>" % k, v).replace("<%s>" % k, v)
        if new != txt:
            # ★ newline="" 必须加：Windows 文本模式写回会把 LF 全变 CRLF，
            #   于是一行 URL 修复在 git diff 里显示成整文件重写（实测 116 行）。
            open(fp, "w", encoding="utf-8", newline="").write(new)
            touched.append(rel)
    if not touched:
        return
    git("add", *touched)
    if git("diff", "--cached", "--quiet", check=False).returncode != 0:
        git("commit", "-q", "-m", commit_msg)
        log(step, "、".join(touched) + " 已填充占位符并追加一次提交")


def zenodo_find_doi(repo_url, rounds):
    """轮询 Zenodo 公开 API，找有没有为这个仓库归档的记录。"""
    q = urllib.parse.quote(repo_url)
    for i in range(rounds):
        try:
            req = urllib.request.Request(
                f"{ZENODO_API}?q={q}&size=20",
                headers={"Accept": "application/json",
                         "User-Agent": "reproducible-code-deposit"})
            with urllib.request.urlopen(req, timeout=45) as r:
                data = json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:
            log("6/8", f"查询 Zenodo 失败（第 {i + 1} 轮）：{type(e).__name__}")
            data = None
        if data:
            needle = repo_url.rstrip("/").lower()
            for hit in data.get("hits", {}).get("hits", []):
                blob = json.dumps(hit, ensure_ascii=False).lower()
                if needle in blob and hit.get("doi"):
                    return hit["doi"], hit.get("links", {}).get("self_html", "")
        if i + 1 < rounds:
            print(f"      … 第 {i + 1}/{rounds} 轮未就绪，30 秒后重试", flush=True)
            time.sleep(30)
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-file", default=os.path.join(HERE, "_push", ".ghtoken"))
    ap.add_argument("--repo-name", default=DEFAULT_REPO)
    ap.add_argument("--r-version", default=None)
    ap.add_argument("--license", default="MIT")
    ap.add_argument("--tag", default="v1.0.0")
    ap.add_argument("--no-create", action="store_true")
    ap.add_argument("--allow-placeholder", action="store_true")
    ap.add_argument("--skip-release", action="store_true")
    ap.add_argument("--via-api", action="store_true",
                    help="github.com:443 不可达时，改用 REST API 推送内容与 tag"
                         "（等同调用 push_via_api.py；只依赖 api.github.com）")
    ap.add_argument("--poll", type=int, default=20)
    args = ap.parse_args()

    global _TOKEN
    _TOKEN = os.environ.get("GITHUB_TOKEN", "").strip()
    if not _TOKEN and os.path.exists(args.token_file):
        with open(args.token_file, encoding="utf-8") as fh:
            _TOKEN = fh.read().strip()
    if not _TOKEN:
        raise SystemExit("没有凭据。请用 --token-file 指定，或设置 GITHUB_TOKEN。")
    if not os.path.isdir(os.path.join(REPO_DIR, ".git")):
        raise SystemExit(f"{REPO_DIR} 不是 git 仓库；先跑 make_repo.py 并 git init/commit。")

    # ------------------------------------------------ 1 凭据
    st, me = gh("GET", "/user")
    if st != 200:
        raise SystemExit(f"凭据无效（HTTP {st}）：{redact(me.get('message'))}")
    owner = me["login"]
    log("1/8", f"凭据有效 —— 账号 {owner}")

    # ------------------------------------------------ 2 仓库
    st, repo = gh("GET", f"/repos/{owner}/{args.repo_name}")
    if st == 200:
        log("2/8", f"仓库已存在：{repo['full_name']}"
                  f"（{'私有' if repo.get('private') else '公开'}）")
        if repo.get("private"):
            print("      ! Zenodo 的 GitHub 集成只对**公开**仓库生效；"
                  "若拿不到 DOI，请先公开该仓库。")
    elif st == 404:
        if args.no_create:
            raise SystemExit("仓库不存在，但指定了 --no-create。")
        st2, created = gh("POST", "/user/repos", {
            "name": args.repo_name,
            "description": "Reproducible analysis code for a public-data (NCBI GEO) study "
                           "re-attributing the adenoma-carcinoma expression shift to "
                           "epithelial and stromal compartments. Companion code for a "
                           "manuscript submitted to BMC Genomics.",
            "private": False, "has_issues": True, "has_wiki": False,
            "has_projects": False, "auto_init": False})
        if st2 == 201:
            repo = created
            log("2/8", f"已创建仓库：{repo['full_name']}")
        else:
            raise SystemExit("\n".join([
                f"无法创建仓库（HTTP {st2}）：{redact(created.get('message'))}",
                "",
                "　fine-grained token 默认没有建仓权限。两种办法任选其一：",
                f"　(A) 打开 https://github.com/new 手动建一个名为 {args.repo_name} 的",
                "　　　public 仓库（**不要**勾 Add a README / .gitignore / license），",
                "　　　然后重跑本脚本（会自动跳过建仓）。注意：若 token 是「Only select",
                "　　　repositories」，还要把新仓库加进 token 的授权列表。",
                "　(B) 换一个 classic token（Settings → Developer settings → Tokens",
                "　　　(classic) → 勾 repo），然后重跑本脚本，建仓与推送一次完成。"]))
    else:
        raise SystemExit(f"探测仓库失败（HTTP {st}）：{redact(repo.get('message'))}")

    html_url = f"https://github.com/{owner}/{args.repo_name}"
    clone_url = f"{html_url}.git"

    # ------------------------------------------------ 3 remote + 填仓库地址
    remotes = git("remote").stdout.split()
    if "origin" in remotes:
        git("remote", "set-url", "origin", clone_url)
    else:
        git("remote", "add", "origin", clone_url)
    log("3/8", f"remote origin -> {clone_url}")

    fill_placeholders({"REPO_URL": html_url},
                      "README: fill in repository URL", "3/8")

    # ------------------------------------------------ 4 push（token 只在命令行，不落盘）
    # ★ 本地 tag 先重打：preflight_release.py 会校验「最新 tag 指向 HEAD」，
    #   且 API 推送要把远端 tag 一起移到新提交上。
    git("tag", "-f", args.tag)

    if args.via_api:
        # ★ api.github.com 与 github.com:443 会**各自独立失效**（国内网络常见）。
        #   这条路只用 REST API，完全绕开 git 协议。
        r = subprocess.run(
            [sys.executable, os.path.join(HERE, "push_via_api.py"),
             "--token-file", args.token_file,
             "--slug", f"{owner}/{args.repo_name}",
             "--branch", "main", "--tag", args.tag],
            cwd=HERE, capture_output=True, text=True,
            encoding="utf-8", errors="replace")
        print((r.stdout or "").rstrip())
        if r.returncode != 0:
            raise SystemExit("API 推送失败：\n" + redact(r.stderr or ""))
        log("4/8", "main 已通过 REST API 推送（绕开 github.com）")
        log("5/8", f"tag {args.tag} 已通过 REST API 更新")
    else:
        tok_url = f"https://x-access-token:{_TOKEN}@github.com/{owner}/{args.repo_name}.git"
        p = git("-c", "credential.helper=", "push", tok_url, "HEAD:refs/heads/main",
                check=False)
        if p.returncode != 0:
            err = redact(p.stderr)
            print("      首次推送被拒，尝试强制推送：")
            print("      " + err.strip().replace("\n", "\n      ")[:600])
            git("-c", "credential.helper=", "push", "--force", tok_url,
                "HEAD:refs/heads/main")
        n = git("rev-list", "--count", "HEAD").stdout.strip()
        log("4/8", f"main 已推送（{n} 个提交）")
        gh("PATCH", f"/repos/{owner}/{args.repo_name}", {"default_branch": "main"})

        # -------------------------------------------- 5 tag
        git("-c", "credential.helper=", "push", "-f", tok_url, f"refs/tags/{args.tag}")
        log("5/8", f"tag {args.tag} 已推送")

    doi = None
    if args.skip_release:
        log("6/8", "已跳过 release（--skip-release）")
    else:
        # -------------------------------------------- 6 release
        # ★ 发 release 是这条链上**唯一不可逆**的动作：Zenodo 会立刻把当时的元数据
        #   冻成一个永久、可被引用的 DOI 记录。所以闸门统一交给 preflight_release.py
        #   跑全套（占位符 / 作者元数据 / 稿件合规 / tag 指向），不再只 grep 一个占位符。
        pf = os.path.join(HERE, "preflight_release.py")
        if os.path.exists(pf) and not args.allow_placeholder:
            r = subprocess.run([sys.executable, pf, "--quiet"], cwd=HERE,
                               capture_output=True, text=True,
                               encoding="utf-8", errors="replace")
            if r.returncode != 0:
                raise SystemExit("\n".join([
                    "已中止：发布前闸门未通过（preflight_release.py）。",
                    "",
                    (r.stdout or r.stderr or "").rstrip(),
                    "",
                    "　清掉上面所有 ⛔ 项再重跑；确实要强行继续，加 --allow-placeholder。",
                    "　（push 与 tag 已完成，重新运行会从 release 这步继续。）"]))
            log("6/8", "发布前闸门已通过")
        elif args.allow_placeholder:
            log("6/8", "警告：--allow-placeholder 已跳过发布前闸门")
        notes = os.path.join(REPO_DIR, "RELEASE_NOTES_v1.0.0.md")
        body = open(notes, encoding="utf-8").read() if os.path.exists(notes) else ""
        st, rel = gh("POST", f"/repos/{owner}/{args.repo_name}/releases", {
            "tag_name": args.tag, "name": args.tag, "body": body,
            "draft": False, "prerelease": False})
        if st == 201:
            log("6/8", f"Release 已发布：{rel['html_url']}")
        elif st == 422 and "already_exists" in json.dumps(rel):
            log("6/8", "Release 已存在，跳过创建")
        else:
            log("6/8", f"Release 创建失败（HTTP {st}）：{redact(rel.get('message'))}")

        # -------------------------------------------- 7 DOI
        if args.poll > 0:
            doi, rec = zenodo_find_doi(html_url, args.poll)
        if doi:
            log("7/8", f"Zenodo DOI 已取回：{doi}")
        else:
            log("7/8", "\n".join([
                "还没等到 Zenodo DOI。常见原因（按概率排序）：",
                "　1. Zenodo 侧还没打开这个仓库的开关 —— 打开",
                "　　 https://zenodo.org/account/settings/github/ 找到该仓库并把开关拨到 ON。",
                "　　 **开关必须在发 release 之前打开**，否则这一次 release 不会被归档；",
                "　　 打开开关后重跑本脚本（会新建一个新版本的 release）即可。",
                "　2. 仓库还是私有 —— Zenodo 只归档公开仓库。",
                "　3. 刚发完 release，Zenodo 还在排队（通常几分钟）。",
                "　稍后单独取 DOI：python push_github.py --skip-release --no-create "
                f"--poll 20"]))

    # ------------------------------------------------ 8 回填稿件
    ca = {"repo_url": html_url, "license": args.license}
    if doi:
        ca["doi"] = f"https://doi.org/{doi}"
    if args.r_version:
        ca["r_version"] = args.r_version
    capath = os.path.join(HERE, "code_availability.json")
    old = {}
    if os.path.exists(capath):
        try:
            old = json.load(open(capath, encoding="utf-8"))
        except Exception:
            old = {}
    old.update(ca)
    with open(capath, "w", encoding="utf-8") as fh:
        json.dump(old, fh, ensure_ascii=False, indent=2)
    log("8/8", f"已写入 code_availability.json：{json.dumps(ca, ensure_ascii=False)}")

    # DOI 占位符一并回填（make_repo.py 打包时也做一次，这一步是拿到 DOI 之后的补充）
    val = ("https://doi.org/" + doi) if doi else \
          "not yet minted (Zenodo assigns one to each tagged release)"
    fill_placeholders({"ARCHIVED_DOI": val},
                      "README: fill in the archived-version DOI", "8/8")

    rc = subprocess.run([sys.executable, os.path.join(HERE, "build_bmc_manuscript.py")],
                        cwd=HERE, capture_output=True, text=True)
    if rc.returncode == 0:
        log("8/8", "投稿稿已重建（Code availability 已用真实值）")
    else:
        log("8/8", "重建稿件失败，请手动运行 build_bmc_manuscript.py：\n"
                   + redact(rc.stderr)[-500:])

    print("\n" + "=" * 72)
    print(f"仓库：{html_url}")
    print(f"DOI ：{('https://doi.org/' + doi) if doi else '（尚未取到）'}")
    if not doi:
        print("      打开 https://zenodo.org/account/settings/github/ 开启该仓库的开关后，")
        print("      重跑：python push_github.py --skip-release --no-create")
    print("=" * 72)


if __name__ == "__main__":
    main()
