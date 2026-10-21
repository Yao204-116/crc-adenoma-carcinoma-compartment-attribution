# -*- coding: utf-8 -*-
"""发布前总闸门：在 Zenodo 拧开关 / 发 GitHub release **之前**跑一遍。

为什么需要它
------------
稿件改错可以再改，仓库推错可以再推；**发 release 是这条链上唯一不可逆的动作** ——
Zenodo 会立刻把当时的元数据冻成一个永久、可被引用的 DOI 记录。所以「能不能发」
不能靠人肉记忆，得写成可执行检查，一次跑全。

它兜住的三类错误（都在实测中真实发生过）
----------------------------------------
1. **占位符漏网**：`<REPO_URL>` 只被填进 README，`RELEASE_NOTES_v1.0.0.md` 里那句
   `git clone <REPO_URL>` 漏了 —— 而 release notes 正是贴进 GitHub release 正文的文字。
   现在改为扫**整个打包产物**，任何文本文件里的占位符都算 FAIL。
2. **内部 TODO 外泄**：`.zenodo.json` 的 `notes` 字段会被原文发布进公开 DOI 记录，
   早期版本写着 "AUTHOR METADATA MUST BE COMPLETED…"，一旦发布就永久挂在那里。
   现在对元数据文本做 TODO 关键词扫描。
3. **作者元数据没填就发**：`creators` 仍是 `__FILL_ME__` → 永久记录里作者是占位符。
   这是最贵的一种，因为 DOI 记录只能通过发「新版本」来修正，旧记录永久存在。

用法
----
    python preflight_release.py            # 全量检查
    python preflight_release.py --quiet    # 只打印失败项

退出码：0 = 可以发；1 = 有阻塞项，不要发。
"""
import argparse
import json
import os
import re
import subprocess
import sys

SRC = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(SRC, "repo")
MS = os.path.join(SRC, "Manuscript_BMC_Genomics_投稿版.html")

# 打包产物里**任何**文本文件都不允许出现的占位符
BAD_TOKENS = ["__FILL_ME__", "<REPO_URL>", "<ARCHIVED_DOI>",
              "\\<REPO_URL\\>", "\\<ARCHIVED_DOI\\>"]
# 元数据文本里不允许出现的内部 TODO 措辞（这些字段会公开）
TODO_WORDS = ["MUST BE COMPLETED", "TODO", "TBD", "待填", "待补", "XXX", "__FILL_ME__"]

# Zenodo legacy deposit schema 的枚举（GitHub 集成读 .zenodo.json 时用它校验）
UPLOAD_TYPES = {"publication", "poster", "presentation", "dataset", "image",
                "video", "software", "lesson", "physicalobject", "other"}
ACCESS_RIGHTS = {"open", "embargoed", "restricted", "closed"}
# 常见合法 license id（legacy vocab）。不在表里只给 WARN，不当 FAIL —— 只列软件/数据常用项。
KNOWN_LICENSES = {"MIT", "GPL-2.0", "GPL-3.0", "Apache-2.0", "BSD-2-Clause",
                  "BSD-3-Clause", "MPL-2.0", "LGPL-2.1", "LGPL-3.0", "ISC",
                  "cc-zero", "cc-by-4.0", "cc-by-sa-4.0", "cc-by-nc-4.0",
                  "other-open", "other-closed"}
RELATION_SCHEMES = {"doi", "handle", "url", "urn", "pmid", "pmcid", "arxiv",
                    "isbn", "issn", "ads", "bibcode", "other"}

# 稿件里的占位符：**不再手写清单**，改成通用规则 —— 任意"内容不以数字开头"的方括号。
# ☆ 为什么改：手写清单抓不到 `[e-mail]` / `[name]` / `[initials]` / `[Initials]`
#   —— 实测旧规则漏了 **8 处**。而正文引用是 `[1]` / `[11, 12]` / `[11-13]`，一律以数字开头，
#   天然被排除；已核 References 区非数字方括号为 **0**，所以通用规则不产生误报。
# ★ 内容上限**不能太小**：DOI 占位符有 118 字符，限 80 会让它整个匹配不到 →
#   报告成"已回填"这种**假通过**（实测踩过）。
# ★ 上限也不能**设死**：{0,300} 意味着"方括号里超过 302 字符就整条漏检"，
#   而这恰恰是一种**最像占位符**的形态 —— 作者把整句说明写进方括号（如 Funding 的
#   "[The authors disclosed receipt of ... identify the grant that supported this work
#   ...]"）。实测 400 字符时命中为 0，闸门会报"已回填"。故放宽到 2000，只留防灾难
#   回溯的保护；真正的边界是 `]`，不是长度。
BRACKET = re.compile(r"\[([^\]\d][^\]]{0,2000})\]")
# DOI 占位符只能等拿到 DOI 之后再填 —— 发布前存在是**预期**的，只提示不阻塞
DOI_SIG = "DOI of the frozen release"

TEXT_EXT = {".md", ".json", ".txt", ".py", ".r", ".csv", ".html", ".tsv",
            ".yml", ".yaml", ".toml", ".cfg", ".ini", ""}

# 这四个脚本的**职责本身就是替换/校验占位符**，源码里出现占位符字面量是设计使然
# （`s.replace("<REPO_URL>", …)`），扫它们只会制造误报。其余文件一律必须干净。
SCAN_EXCLUDE = {"make_repo.py", "push_github.py", "push_via_api.py",
                "preflight_release.py", "fill_author_metadata.py"}

results = []          # (severity, name, ok, detail)


def add(sev, name, ok, detail=""):
    results.append((sev, name, ok, detail))


def block(name, ok, detail=""):
    add("BLOCK", name, ok, detail)


def warn(name, ok, detail=""):
    add("WARN", name, ok, detail)


# ---------------------------------------------------------------- A 元数据闸门
meta_path = os.path.join(REPO, ".zenodo.json")
meta = None
if not os.path.exists(meta_path):
    block("`.zenodo.json` 存在", False, meta_path)
else:
    try:
        meta = json.load(open(meta_path, encoding="utf-8"))
        block("`.zenodo.json` 可解析", True)
    except Exception as e:
        block("`.zenodo.json` 可解析", False, str(e))

if meta is not None:
    # A1 必填字段
    need = ["title", "upload_type", "access_right", "license", "description",
            "creators"]
    lack = [k for k in need if not str(meta.get(k, "")).strip()]
    block("元数据必填字段齐全", not lack, f"缺 {lack}" if lack else "")

    # A2 ★ 作者元数据（最贵的一项：占位符会被冻进永久 DOI 记录）
    cr = meta.get("creators") or []
    bad = []
    for i, c in enumerate(cr):
        if not isinstance(c, dict):
            bad.append(f"#{i} 不是对象")
            continue
        nm = str(c.get("name", ""))
        af = str(c.get("affiliation", ""))
        if (not nm.strip()) or "__FILL_ME__" in nm:
            bad.append(f"#{i} name={nm!r}")
        if (not af.strip()) or "__FILL_ME__" in af:
            bad.append(f"#{i} affiliation={af!r}")
        oid = str(c.get("orcid", "")).strip()
        if oid and not re.match(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$", oid):
            bad.append(f"#{i} orcid 格式={oid!r}")
    block("作者元数据已填（creators 非占位符）", bool(cr) and not bad,
          ("、".join(bad) if bad else f"{len(cr)} 位作者"))

    # A3 枚举
    ut = meta.get("upload_type", "")
    block("upload_type 合法", ut in UPLOAD_TYPES, f"{ut!r}")
    ar = meta.get("access_right", "")
    block("access_right 合法", ar in ACCESS_RIGHTS, f"{ar!r}")
    lic = meta.get("license", "")
    block("license 合法", lic in KNOWN_LICENSES, f"{lic!r}")
    lang = meta.get("language", "")
    warn("language 为 ISO 639-3 三字母码", len(lang) == 3, f"{lang!r}")

    # A4 关联标识符：scheme 必须显式（不写默认按 doi 解析 → 自由文本会被判非法）
    ri = meta.get("related_identifiers") or []
    ri_bad = []
    for i, r in enumerate(ri):
        if not isinstance(r, dict):
            ri_bad.append(f"#{i} 不是对象")
            continue
        sch = r.get("scheme")
        ident = str(r.get("identifier", "")).strip()
        if sch is None:
            ri_bad.append(f"#{i} 缺 scheme（会默认按 doi 校验）")
        elif sch not in RELATION_SCHEMES:
            ri_bad.append(f"#{i} scheme={sch!r} 不在枚举内")
        if sch == "doi" and re.match(r"^https?://", ident):
            ri_bad.append(f"#{i} doi 域须写裸 DOI，不带 https:// 前缀")
        if not r.get("relation"):
            ri_bad.append(f"#{i} 缺 relation")
    block("related_identifiers 的 scheme 显式且合法", not ri_bad,
          "；".join(ri_bad) if ri_bad else f"{len(ri)} 条")

    # A5 ★ 元数据文本不得含内部 TODO（这些字段会原文公开）
    blob = " ".join(str(meta.get(k, "")) for k in
                    ("title", "description", "notes", "version"))
    blob += " " + " ".join(map(str, meta.get("keywords", [])))
    blob += " " + " ".join(str(c.get("name", "")) + str(c.get("affiliation", ""))
                           for c in cr if isinstance(c, dict))
    hits = [w for w in TODO_WORDS if w.lower() in blob.lower()]
    block("元数据文本无内部 TODO 措辞", not hits, f"命中 {hits}" if hits else "")

# ---------------------------------------------------------------- B 稿件闸门
if not os.path.exists(MS):
    block("BMC 投稿版稿件存在", False, os.path.basename(MS))
    ms_html = ""
else:
    block("BMC 投稿版稿件存在", True)
    ms_html = open(MS, encoding="utf-8").read()

if ms_html:
    ms_text = re.sub(r"<[^>]+>", " ", re.sub(r"<style>.*?</style>", "", ms_html,
                                             flags=re.S))

    def _tokens(text):
        """取出非数字方括号占位符的**完整**内容 `[…]`。

        两条经验：
        * 不能按固定字数切片 —— 原文里占位符后面紧跟换行与上标数字，切 70 字会把
          下一段文字一起吞进来，报告没法读（实测如此）；
        * 不能手写关键词清单 —— 清单永远漏（见 BRACKET 处的说明）。"""
        return [re.sub(r"\s+", " ", m.group(0)).strip()
                for m in BRACKET.finditer(text)]

    def _fmt(lst):
        seen = {}
        for t in lst:
            seen[t] = seen.get(t, 0) + 1
        return "；".join(f"{k} ×{v}" if v > 1 else k for k, v in seen.items())

    _toks = _tokens(ms_text)
    left_b = [t for t in _toks if DOI_SIG not in t]
    left_p = [t for t in _toks if DOI_SIG in t]
    block("稿件无任何未填占位符（任意非数字方括号）", not left_b,
          (f"{len(left_b)} 处：" + _fmt(left_b)) if left_b else "")
    warn("稿件 DOI 占位符（拿到 DOI 后回填，发布前允许存在）",
         not left_p, f"{len(left_p)} 处待回填" if left_p else "已回填")

    # B3 借用 BMC 自检做整体合规（21 项）
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    try:
        p = subprocess.run([sys.executable, "verify_bmc.py"], cwd=SRC, env=env,
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=180)
        out = (p.stdout or "") + (p.stderr or "")
        n_pass = len(re.findall(r"^\[PASS\]", out, re.M))
        n_fail = len(re.findall(r"^\[FAIL\]", out, re.M))
        fails = re.findall(r"^\[FAIL\] (.*)$", out, re.M)
        block("BMC 自检全部通过", n_fail == 0 and "全部通过" in out,
              f"{n_pass} 通过 / {n_fail} 失败" +
              ("：" + "；".join(fails) if fails else ""))
    except Exception as e:
        block("BMC 自检全部通过", False, f"无法运行 verify_bmc.py: {e}")

# ---------------------------------------------------------------- C 打包产物闸门
if not os.path.isdir(REPO):
    block("打包产物 repo/ 存在", False, REPO)
else:
    block("打包产物 repo/ 存在", True)
    hits = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d != ".git"]
        for f in files:
            ext = os.path.splitext(f)[1].lower()
            if ext not in TEXT_EXT or f in SCAN_EXCLUDE:
                continue
            fp = os.path.join(root, f)
            try:
                s = open(fp, encoding="utf-8", errors="ignore").read()
            except Exception:
                continue
            for tok in BAD_TOKENS:
                if tok in s:
                    rel = os.path.relpath(fp, REPO).replace("\\", "/")
                    for i, line in enumerate(s.splitlines(), 1):
                        if tok in line:
                            hits.append(f"{rel}:{i} `{tok}`")
    block("★ 打包产物无任何占位符残留", not hits,
          f"{len(hits)} 处：" + "；".join(hits[:6]) if hits else "")

    # C2 代码可用性回填
    ca_path = os.path.join(SRC, "code_availability.json")
    ca = {}
    if os.path.exists(ca_path):
        try:
            ca = json.load(open(ca_path, encoding="utf-8"))
        except Exception:
            ca = {}
    lack = [k for k in ("repo_url", "license", "r_version") if not ca.get(k)]
    block("code_availability.json 三键齐全", not lack, f"缺 {lack}" if lack else "")
    warn("code_availability.json 已含 doi（发布后回填）", bool(ca.get("doi")),
         ca.get("doi") or "尚无 DOI —— 取回后回填")

    # C3 版本号与 tag 一致
    if meta is not None:
        ver = str(meta.get("version", "")).strip()
        try:
            tags = subprocess.run(["git", "-C", REPO, "tag"], capture_output=True,
                                  text=True).stdout.split()
        except Exception:
            tags = []
        ok = ver in tags
        block("`.zenodo.json` 的 version 与某个 tag 一致", ok,
              f"version={ver!r}；tags={tags}")

    # ------------------------------------------------------------ D git 闸门
    def _git(*a):
        try:
            return subprocess.run(["git", "-C", REPO] + list(a), capture_output=True,
                                  text=True).stdout.strip()
        except Exception:
            return ""

    dirty = _git("status", "--porcelain")
    block("git 工作区干净（无未提交改动）", not dirty,
          dirty.splitlines()[0] if dirty else "")
    tags = _git("tag").split()
    block("至少存在一个 tag", bool(tags), f"tags={tags}")
    head = _git("rev-parse", "HEAD")
    if tags:
        # ★ 真实的约束是「**最新那个** tag 指向 HEAD」，历史 tag 停在旧提交是**正常**的
        #   （发过 v1.0.0 后升到 v1.0.1，v1.0.0 当然不等于 HEAD）。
        #   旧实现逐个要求全部 tag == HEAD，一旦发第二个版本就永远拦不通过。
        #   "最新"必须是**语义版本序**（1.0.10 > 1.0.9），不能按字典序 ——
        #   字典序会把 v1.0.1 排在 v1.0.0 **前面**，于是报出「v1.0.0->xxx(≠HEAD)」。
        def _vkey(t):
            m = re.match(r"^v?(\d+(?:\.\d+)*)", t)
            if not m:
                return ((), -1, t)      # 非版本号形态的 tag 排最后，仅用于稳定输出
            nums = tuple(int(x) for x in m.group(1).split("."))
            return (nums, 0, t)

        try:
            newest = max(tags, key=_vkey)
        except Exception:
            newest = tags[-1]
        tcs = {t: _git("rev-list", "-n", "1", t) for t in tags}
        info = [f"{t}->{tcs[t][:7]}{'' if tcs[t] == head else '(≠HEAD)'}"
                for t in sorted(tags, key=_vkey)]
        # ★ 发布后的正常状态：DOI 回填必然产生「tag 之后的提交」（README 加 DOI 行、
        #   .zenodo.json 同步）。此时 tag 是 HEAD 的**祖先**，内容没变，不该拦。
        #   旧实现只认 tag == HEAD → 每次回填后闸门永久 FAIL；一个总在报警的闸门等于没有闸门。
        #   判据：HEAD..tag 的提交数为 0 ⟺ tag 已被 HEAD 完整包含。
        newest_sha = tcs.get(newest, "")
        behind = _git("rev-list", "--count", "HEAD..%s" % newest) if newest else ""
        ahead = _git("rev-list", "--count", "%s..HEAD" % newest) if newest else ""
        is_ancestor = (behind == "0")
        detail = "最新=%s ｜ %s" % (newest, " ｜ ".join(info))
        if newest_sha == head:
            pass
        elif is_ancestor:
            detail += " ｜ HEAD 比 tag 多 %s 个提交（发布后回填，正常）" % (ahead or "?")
        block("最新 tag 已包含在 HEAD 中（== 或为其祖先）",
              newest_sha == head or is_ancestor, detail)

# ---------------------------------------------------------------- 报告
ap = argparse.ArgumentParser(description="发布前总闸门")
ap.add_argument("--quiet", action="store_true", help="只打印未通过项（供脚本内调用）")
ARGS = ap.parse_args()

nb = sum(1 for s, _n, ok, _d in results if s == "BLOCK" and not ok)
nw = sum(1 for s, _n, ok, _d in results if s == "WARN" and not ok)
print("=" * 74)
print("发布前闸门 · preflight_release")
print("=" * 74)
for sev, name, ok, detail in results:
    if ok:
        if ARGS.quiet:
            continue
        mark, tag = "PASS", ""
    else:
        mark = "FAIL" if sev == "BLOCK" else "WARN"
        tag = "  ⛔ 阻塞" if sev == "BLOCK" else "  ⚠ 提示"
    line = f"[{mark}] {name}"
    if detail:
        line += f"  — {detail}"
    print(line + tag)
print("-" * 74)
print(f"阻塞项 {nb} 个 ｜ 提示项 {nw} 个")
if nb:
    print("\n→ 不要发 release。先清掉上面所有 ⛔ 项。")
    print("   （作者/邮箱/基金/CRediT 只能由本人提供；DOI 只能等 release 之后回填。）")
else:
    print("\n→ 可以发 release。发完记得：")
    print("   1) 回填 code_availability.json 的 doi 与稿件 DOI 占位符；")
    print("   2) 重跑 make_repo.py；本闸门会转为提示「DOI 已回填」。")
print("=" * 74)
sys.exit(1 if nb else 0)
