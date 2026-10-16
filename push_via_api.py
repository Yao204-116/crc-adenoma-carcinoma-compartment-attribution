# -*- coding: utf-8 -*-
"""当 `github.com:443` 被阻断（而 `api.github.com` 仍可达）时，用 REST API 推送。

**为什么需要**：`git push` 走 `github.com:443`；REST API 走 `api.github.com`。
国内网络下这两个域名经常一个通、一个不通 —— 实测同一分钟内
`api.github.com` 返回 200（0.38 s），而 `github.com` 连接超时（21 s）、
`git ls-remote` 直接失败。此时 git 协议完全用不了，但 API 可以写。

**做法**：不用 git 协议，改用 GitHub Git Data API 构造一个与本地 HEAD 等价的提交：
  1. 取远端分支当前 commit 及其**完整树**（recursive）→ path → blob sha
  2. 取本地树（`git ls-tree -r -z HEAD`）→ 逐 path 比较，**只对差异文件**建 blob
  3. `POST /git/trees` 用 `base_tree` + 差异条目（新增/修改给 `content`，删除给 `sha: null`）
  4. `POST /git/commits`（message 取本地 HEAD 的 message，parent = 远端当前 commit）
  5. `PATCH /git/refs/heads/<branch>` 快进分支；`PATCH /git/refs/tags/<tag>` 强制更新 tag

只上传**差异**，所以即使仓库有 15 MB 也只有几百字节的请求体。

用法::

    python push_via_api.py --token-file _push/.ghtoken              # 同步 main + tag
    python push_via_api.py --token-file _push/.ghtoken --dry-run    # 只看差异
    python push_via_api.py --token-file _push/.ghtoken --no-tag

★ 凭据只从文件读、只放进 Authorization 头，**绝不打印**；所有回显过 `redact()`。
"""
import argparse
import base64
import json
import os
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.join(HERE, "repo")
API = "https://api.github.com"
_TOK = None


def redact(s):
    """把 token 从任何要打印的文本里抹掉。"""
    if not s:
        return s
    out = str(s)
    if _TOK:
        out = out.replace(_TOK, "***REDACTED***")
        out = out.replace(urllib.parse.quote(_TOK, safe=""), "***REDACTED***")
    return out


def api(method, path, payload=None):
    """调 GitHub REST API；返回 (status, 解析后的 JSON 或原始文本)。"""
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    req = urllib.request.Request(API + path, data=data, method=method)
    req.add_header("Authorization", "Bearer %s" % _TOK)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("X-GitHub-Api-Version", "2022-11-28")
    req.add_header("User-Agent", "reproducible-code-deposit")
    if data is not None:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=90) as r:
            raw = r.read().decode("utf-8", "replace")
            return r.status, (json.loads(raw) if raw else {})
    except urllib.error.HTTPError as e:
        raw = e.read().decode("utf-8", "replace")
        try:
            return e.code, json.loads(raw)
        except Exception:
            return e.code, raw


def git(*args):
    return subprocess.run(["git"] + list(args), cwd=REPO_DIR,
                          capture_output=True, text=True)


def local_tree():
    """本地 HEAD 的 path → blob sha。"""
    p = git("ls-tree", "-r", "-z", "HEAD")
    if p.returncode != 0:
        raise SystemExit("git ls-tree 失败：\n" + p.stderr)
    out = {}
    for rec in p.stdout.split("\0"):
        if not rec:
            continue
        meta, path = rec.split("\t", 1)
        mode, typ, sha = meta.split()
        if typ == "blob":
            out[path] = sha
    return out


def local_blob(path):
    """按 git 的存储形态取出文件内容（走 git cat-file，自动处理换行归一化）。"""
    p = subprocess.run(["git", "show", "HEAD:" + path], cwd=REPO_DIR,
                       capture_output=True)
    if p.returncode != 0:
        raise SystemExit("git show 失败：" + path + "\n" + p.stderr.decode("utf-8", "replace"))
    return p.stdout


def local_commit_meta():
    """本地 HEAD 的作者/提交者身份与时间。

    ★ 必须原样带给 API：这样 GitHub 造出的 commit 与本地**逐字段相同**，
      sha 也就会相同 —— 否则远端会多出一个内容一样但 sha 不同的提交，
      下次 git 网络恢复后 `git push` 会被判 non-fast-forward 而分叉。
    """
    f = "%an%x00%ae%x00%aI%x00%cn%x00%ce%x00%cI"
    p = git("log", "-1", "--pretty=" + f)
    if p.returncode != 0:
        raise SystemExit("读本地提交元数据失败：\n" + p.stderr)
    an, ae, ad, cn, ce, cd = p.stdout.strip().split("\x00")
    return ({"name": an, "email": ae, "date": ad},
            {"name": cn, "email": ce, "date": cd})


def remote_tree(sha):
    """远端某 commit 的树 → path → blob sha（recursive）。"""
    st, d = api("GET", "/repos/%s/git/trees/%s?recursive=1" % (SLUG, sha))
    if st != 200:
        raise SystemExit("取远端树失败 (%s)：%s" % (st, d))
    return {e["path"]: e["sha"] for e in d.get("tree", []) if e["type"] == "blob"}


SLUG = None
BRANCH = None


def sync_tag(tag, sha):
    """把远端 tag 指到 sha（仅支持轻量 tag）。

    Zenodo 归档的是 **tag 指向的那个 commit**，所以打完 tag 又改动文件时
    必须重打 —— 否则 DOI 记录对应的是旧内容。
    """
    st, tg = api("GET", "/repos/%s/git/ref/tags/%s" % (SLUG, tag))
    if st == 200 and tg["object"]["type"] == "commit":
        st2, _r2 = api("PATCH", "/repos/%s/git/refs/tags/%s" % (SLUG, tag),
                       {"sha": sha, "force": True})
        print("      tag %s %s → %s" % (tag,
              "已更新" if st2 in (200, 201) else "更新失败(%s)" % st2, sha[:8]))
    elif st == 200:
        print("      tag %s 是 annotated tag（非轻量），本脚本不重建 —— "
              "请在网络恢复后 `git tag -f` 重推。" % tag)
    else:
        print("      tag %s 不存在于远端，跳过。" % tag)


def main():
    global _TOK, SLUG, BRANCH
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-file", default=os.path.join(HERE, "_push", ".ghtoken"))
    ap.add_argument("--slug", default="Yao204-116/crc-adenoma-carcinoma-compartment-attribution")
    ap.add_argument("--branch", default=None)
    ap.add_argument("--tag", default="v1.0.0")
    ap.add_argument("--no-tag", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    _TOK = open(args.token_file, encoding="utf-8").read().strip()
    SLUG = args.slug

    if not os.path.isdir(os.path.join(REPO_DIR, ".git")):
        raise SystemExit("%s 不是 git 仓库" % REPO_DIR)

    # 本地 HEAD
    head = git("rev-parse", "HEAD").stdout.strip()
    msg = git("log", "-1", "--pretty=%B").stdout.rstrip("\n")
    br = args.branch or git("rev-parse", "--abbrev-ref", "HEAD").stdout.strip()
    BRANCH = br

    st, repo = api("GET", "/repos/%s" % SLUG)
    print("[1/5] 仓库 %s（默认分支 %s）" % (repo.get("full_name"), repo.get("default_branch")))

    st, ref = api("GET", "/repos/%s/git/ref/heads/%s" % (SLUG, br))
    if st != 200:
        raise SystemExit("取远端分支失败 (%s)：%s" % (st, ref))
    remote_head = ref["object"]["sha"]
    print("[2/5] 远端 %s = %s ；本地 HEAD = %s" % (br, remote_head[:8], head[:8]))
    if remote_head == head:
        # ★ 内容一致时**不能直接 return**：本地可能重打了 tag（`git tag -f`），
        #   此时远端 tag 还指着旧提交。只推内容、不管 tag 会让 release 归档到错的
        #   那个 commit 上（Zenodo 归档的是 tag 指向的 commit）。
        print("      内容已一致，跳过内容推送（仍同步 tag）。")
        if not args.no_tag and not args.dry_run:
            sync_tag(args.tag, head)
        if args.dry_run:
            print("\n--dry-run：未做任何写操作。")
        return

    rt = remote_tree(remote_head)
    lt = local_tree()
    added = sorted(set(lt) - set(rt))
    removed = sorted(set(rt) - set(lt))
    changed = sorted(p for p in set(lt) & set(rt) if lt[p] != rt[p])
    print("[3/5] 差异：新增 %d  修改 %d  删除 %d" % (len(added), len(changed), len(removed)))
    for p in added:
        print("      + %s" % p)
    for p in changed:
        print("      M %s" % p)
    for p in removed:
        print("      - %s" % p)
    if args.dry_run:
        print("\n--dry-run：未做任何写操作。")
        return

    entries = []
    for p in added + changed:
        entries.append({"path": p, "mode": "100644", "type": "blob",
                        "content": local_blob(p).decode("utf-8", "surrogateescape")})
    for p in removed:
        entries.append({"path": p, "mode": "100644", "type": "blob", "sha": None})

    author, committer = local_commit_meta()

    st, rc = api("GET", "/repos/%s/git/commits/%s" % (SLUG, remote_head))
    if st != 200:
        raise SystemExit("取远端 commit 失败 (%s)：%s" % (st, rc))
    base_tree = rc["tree"]["sha"]

    st, tree = api("POST", "/repos/%s/git/trees" % SLUG,
                   {"base_tree": base_tree, "tree": entries})
    if st not in (200, 201):
        raise SystemExit("建树失败 (%s)：%s" % (st, tree))
    print("[4/5] 新树 %s（base_tree %s）" % (tree["sha"][:8], base_tree[:8]))

    st, commit = api("POST", "/repos/%s/git/commits" % SLUG,
                     {"message": msg, "tree": tree["sha"], "parents": [remote_head],
                      "author": author, "committer": committer})
    if st not in (200, 201):
        raise SystemExit("建提交失败 (%s)：%s" % (st, commit))
    same = commit["sha"] == head
    print("      新提交 %s%s" % (commit["sha"][:8],
          "  ✓ 与本地 HEAD 逐字段一致（sha 相同，不会分叉）" if same
          else "  ⚠ 与本地 HEAD 不同（树条目顺序等差异），本地需 fetch 对齐"))

    st, upd = api("PATCH", "/repos/%s/git/refs/heads/%s" % (SLUG, br),
                  {"sha": commit["sha"], "force": True})
    if st not in (200, 201):
        raise SystemExit("更新分支失败 (%s)：%s" % (st, upd))
    print("[5/5] %s 已更新 → %s" % (br, commit["sha"][:8]))

    if not args.no_tag:
        sync_tag(args.tag, commit["sha"])

    print("\n完成。注意：这一步绕开了 git 协议，**远端提交的 sha 与本地不同** —— "
          "下次网络恢复后建议先 `git fetch` 再对齐。")


if __name__ == "__main__":
    main()
