# -*- coding: utf-8 -*-
"""把作者元数据从一个 JSON 分发到「仓库 .zenodo.json」与「稿件标题页 / Declarations」。

为什么需要这个脚本：

  发 release 是这条链上**唯一不可逆**的动作 —— Zenodo 会立刻把当时的元数据冻成一个
  永久、可被引用的 DOI 记录。而作者信息天然出现在两个地方：

    1. 仓库根 `.zenodo.json` 的 creators（真源在 `_repo_assets/.zenodo.json`）
    2. 稿件标题页 / Funding / Authors' contributions

  手填两处必然漂移（改了一处忘了另一处），所以只维护 `author_metadata.json` 一个文件，
  由本脚本分发。`build_bmc_manuscript.py` 会读同一个文件生成标题页与 Declarations。

用法：
    python fill_author_metadata.py --init     # 生成可编辑的模板 author_metadata.json
    python fill_author_metadata.py --check    # 只校验，不写任何东西
    python fill_author_metadata.py            # 校验 → 写 .zenodo.json → 重建稿件 → 跑闸门
"""
import argparse
import json
import os
import re
import subprocess
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)

META_PATH = "author_metadata.json"
ZENODO_SRC = os.path.join("_repo_assets", ".zenodo.json")
ZENODO_DST = os.path.join("repo", ".zenodo.json")

ORCID_RE = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$")

TEMPLATE = {
    "_howto": [
        "authors 按稿件署名顺序排列。aff 是单位编号，对应 affiliations 的键。",
        "corresponding 只能有一个为 true；该作者的邮箱写在顶层 corresponding.email。",
        "orcid 可以不填；填的话用裸编号（0000-0002-1825-0097），不要带 https://orcid.org/ 前缀。",
        "initials 用于 Authors' contributions（如 Yao Zhang -> YZ）。credit 为空时不会用到。",
        "credit 的键是 initials，值是 CRediT 角色名列表；不确定就留 {}，脚本会退回占位提示。",
        "funding 写整句英文。无基金就写：This research did not receive any specific "
        "grant from funding agencies in the public, commercial, or not-for-profit sectors.",
        "★ 填完先跑 --check，通过后再跑正式模式（会重建稿件并写 .zenodo.json）。",
    ],
    "authors": [
        {"name": "", "aff": [1], "corresponding": True, "orcid": "", "initials": ""},
    ],
    "affiliations": {"1": ""},
    "corresponding": {"email": ""},
    "funding": "",
    "credit": {},
}


def die(msg):
    print("✗ " + msg)
    sys.exit(1)


def load():
    if not os.path.exists(META_PATH):
        die("%s 不存在。先跑 --init 生成模板。" % META_PATH)
    try:
        return json.load(open(META_PATH, encoding="utf-8"))
    except Exception as e:
        die("%s 不是合法 JSON：%s" % (META_PATH, e))


def validate(m):
    """返回 (hard, soft)。

    ★ 区分「填错了」与「还没填」，因为两者的正确处置完全不同：
      * hard（格式错误）→ **任何模式都拒绝写盘**。空姓名/坏 ORCID/引用不存在的单位编号
        属于这一类：写进去只会产出比占位符更糟的东西（空名字、坏元数据）。
      * soft（只是还没填）→ 默认也拒绝（宁缺勿滥），但 `--partial` 允许先写已完成的部分。
        邮箱缺失会渲染成可见的 `[e-mail]` 占位，而闸门会把它抓出来 —— 所以不会误发。
    """
    hard, soft = [], []
    authors = m.get("authors")
    if not isinstance(authors, list) or not authors:
        hard.append("authors 必须是非空列表")
        return hard, soft
    affs = m.get("affiliations") or {}
    if not affs:
        hard.append("affiliations 不能为空")

    n_corr = 0
    for i, a in enumerate(authors):
        tag = "authors[%d]" % i
        nm = (a.get("name") or "").strip()
        if not nm:
            hard.append("%s.name 为空" % tag)
        if "_FILL_ME_" in nm or "__FILL_ME__" in nm:
            hard.append("%s.name 仍是占位符" % tag)
        a_aff = a.get("aff") or []
        if not a_aff:
            hard.append("%s.aff 为空" % tag)
        for k in a_aff:
            if str(k) not in {str(x) for x in affs}:
                hard.append("%s.aff 引用了不存在的单位编号 %s" % (tag, k))
        o = (a.get("orcid") or "").strip()
        if o:
            if o.startswith("http"):
                hard.append("%s.orcid 请用裸编号，不要带 URL 前缀" % tag)
            elif not ORCID_RE.match(o):
                hard.append("%s.orcid 格式不对（应为 0000-0002-1825-0097 形式）：%s" % (tag, o))
        if a.get("corresponding"):
            n_corr += 1
        if not (a.get("initials") or "").strip():
            soft.append("%s.initials 为空（credit 非空时会用到）" % tag)

    if n_corr != 1:
        hard.append("corresponding=true 的作者必须恰好 1 个，当前 %d 个" % n_corr)

    corr = m.get("corresponding") or {}
    if "@" not in (corr.get("email") or ""):
        soft.append("corresponding.email 缺失或不含 @ → 稿件会渲染成 [e-mail] 占位")

    for k, v in affs.items():
        if not (v or "").strip():
            hard.append("affiliations[%s] 为空" % k)

    if not (m.get("funding") or "").strip():
        soft.append("funding 为空 → 稿件 Funding 会退回占位提示")

    credit = m.get("credit") or {}
    inis = {(a.get("initials") or "").strip() for a in authors}
    for ini in credit:
        if ini not in inis:
            hard.append("credit 里的 %s 与任何作者的 initials 都不匹配" % ini)

    return hard, soft


def creators_of(m):
    affs = m["affiliations"]
    out = []
    for a in m["authors"]:
        aff = "; ".join(str(affs[str(k)]) for k in (a.get("aff") or []))
        c = {"name": a["name"].strip(), "affiliation": aff}
        o = (a.get("orcid") or "").strip()
        if o:
            c["orcid"] = o
        out.append(c)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", action="store_true", help="生成模板 author_metadata.json")
    ap.add_argument("--check", action="store_true", help="只校验，不写任何文件")
    ap.add_argument("--no-rebuild", action="store_true", help="只写 .zenodo.json，不重建稿件")
    ap.add_argument("--no-commit", action="store_true",
                    help="不自动 git commit / 重打 tag（默认会做，否则闸门必拦）")
    ap.add_argument("--partial", action="store_true",
                    help="允许在仍有字段未填时写入已完成的部分（release 仍会被闸门拦住）")
    args = ap.parse_args()

    if args.init:
        if os.path.exists(META_PATH):
            die("%s 已存在，未覆盖。要重来请先自己备份并删除。" % META_PATH)
        with open(META_PATH, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(TEMPLATE, fh, ensure_ascii=False, indent=2)
            fh.write("\n")
        print("已生成模板：%s" % META_PATH)
        for line in TEMPLATE["_howto"]:
            print("  · " + line)
        print("\n填完后：python fill_author_metadata.py --check")
        return

    m = load()
    hard, soft = validate(m)

    print("─" * 66)
    for e in hard:
        print("✗ %s" % e)
    for s in soft:
        print("⚠ 未填：%s" % s)
    if not hard and not soft:
        print("✓ 全部字段齐备")
    print("─" * 66)

    if args.check:
        print("格式错误 %d ｜ 未填 %d" % (len(hard), len(soft)))
        print("--check：未写任何文件。")
        return
    if hard:
        die("%d 处**格式错误**（不是「还没填」，是「填错了」）—— 未写任何文件。"
            % len(hard))
    if soft and not args.partial:
        die("%d 项尚未填写。补齐再跑；或用 --partial 只写入已完成的部分"
            "（草稿会带可见占位符，闸门仍会拦住 release，不会误发）。" % len(soft))

    names = ", ".join(a["name"].strip() for a in m["authors"])
    print("✓ 校验通过：%d 位作者 —— %s" % (len(m["authors"]), names))
    for a in m["authors"]:
        flag = "（通讯）" if a.get("corresponding") else ""
        print("    %-28s 单位 %s %s" % (a["name"].strip(), a.get("aff"), flag))
    if soft:
        print("  ⚠ --partial：本次只写入已完成部分，下列仍是占位 → **发不了 release**：")
        for s in soft:
            print("     · %s" % s)

    # ---------------------------------------------------------- 1 写 .zenodo.json（真源）
    if not os.path.exists(ZENODO_SRC):
        die("找不到 %s" % ZENODO_SRC)
    z = json.load(open(ZENODO_SRC, encoding="utf-8"))
    z["creators"] = creators_of(m)
    with open(ZENODO_SRC, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(z, fh, ensure_ascii=False, indent=2)
        fh.write("\n")
    print("\n[1/4] 已写入 %s（creators %d 位）" % (ZENODO_SRC, len(z["creators"])))

    if args.no_rebuild:
        print("（--no-rebuild：未重建稿件、未跑闸门。记得手动跑 make_repo.py。）")
        return

    # ---------------------------------------------------------- 2 重建稿件 + 重打包
    for script, what in (("build_bmc_manuscript.py", "投稿稿"),
                         ("make_repo.py", "仓库产物")):
        r = subprocess.run([sys.executable, script], capture_output=True, text=True,
                           encoding="utf-8", errors="replace")
        ok = r.returncode == 0
        print("[2/4] %-26s %s" % (script, "OK" if ok else "失败"))
        if not ok:
            print((r.stdout or "")[-800:] + (r.stderr or "")[-800:])
            die("%s 运行失败，先修好再继续。" % script)

    # ---------------------------------------------------------- 3 提交并重打 tag
    # ★ 闸门要求「工作区干净」+「最新 tag 指向 HEAD」。填充必然造成改动，所以这一步
    #   必须一起做掉 —— 否则闸门一定拦（而且拦得对：带着未提交改动的 release 是不该发的）。
    version = z.get("version")
    if not args.no_commit:
        subprocess.run(["git", "-C", "repo", "add", "-A"], capture_output=True)
        r = subprocess.run(["git", "-C", "repo", "commit", "-q", "-m",
                            "Fill author metadata (manuscript title page + .zenodo.json creators)"],
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            print("[3/4] 无改动可提交（或提交失败）：%s" % (r.stdout or r.stderr).strip()[:200])
        head = subprocess.run(["git", "-C", "repo", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
        print("[3/4] 已提交 %s" % head)
        if version:
            subprocess.run(["git", "-C", "repo", "tag", "-f", version],
                           capture_output=True, text=True)
            print("      本地 tag %s 已重打到 %s" % (version, head))
    else:
        print("[3/4] --no-commit：未提交（闸门会因「工作区不干净 / tag 未指向 HEAD」拦截）")

    # ---------------------------------------------------------- 4 跑闸门
    print("[4/4] 发布前闸门 ——")
    r = subprocess.run([sys.executable, "preflight_release.py"], capture_output=True,
                       text=True, encoding="utf-8", errors="replace")
    print((r.stdout or "").rstrip() or (r.stderr or "").rstrip())
    if r.returncode == 0:
        print("\n下一步：发 release（github.com 不通时加 --via-api，会用 REST API 推内容与 tag）：")
        print("    python push_github.py --token-file _push/.ghtoken --via-api")
    else:
        print("\n闸门未通过 —— 上面还有阻塞项。")


if __name__ == "__main__":
    main()
