# -*- coding: utf-8 -*-
"""把复现链打包成可直接推送到 GitHub 的仓库目录。

原则：**只复制，绝不修改原文件**（两个脚本的路径可移植化改造已在源文件完成）。
产物：`repo/`
"""
import json
import os
import shutil
SRC = os.path.dirname(os.path.abspath(__file__))
DST = os.path.join(SRC, "repo")

# ---------------------------------------------------------------- 清单
# 分析脚本（Stage A–I，扁平放置，保持脚本间的相对路径引用有效）
SCRIPTS = """consensus_signature sig_validate drift_diag lcm_diag three_layer_check r3_sens_paneth r3_facts
sc161277_scan sc161277_pairs sc303070_compo sc303070_compo2 sc303070_compo3 sc303070_niche
enrichr_dual gse117606_paired build_suppl_tables make_figs make_figs_sc make_figS1 make_figS3
geo_dl geo_dl2 dl_ebi dl_loop dl_persist sc161277_dl_final gse303070_probe gse303070_tissue
build_bmc_manuscript _fetch_bmc_refs verify_bmc gen_bmc_readiness make_upload_package
make_repo push_github push_via_api preflight_release
fig_en check_fig_english""".split()
# ★ fig_en 是**必须**进仓库的：它是显示名的唯一英文映射，
#   make_figs* / build_suppl_tables 都 import 它；缺了它复现链直接断。
# ★ push_via_api 是 push_github 的兜底：github.com:443 被阻断而 api.github.com 可达时用。
# ★ preflight_release 是**发 release 之前**的总闸门：发布是唯一不可逆的动作，
#   它把「能不能发」写成可执行检查（占位符 / 元数据 / 自检 / tag 指向）。

# 落盘的中间产物（数字来源，供复算核对）
MIDDLE = """panels_def.json suppl_table1_313genes.csv manuscript_facts.json r3_facts.json
r3_sens_paneth.json three_layer_pattern.json three_layer_results.json sc161277_scan.json
sc161277_pairs.json drift_diag.json sc303070_niche.json sc303070_composition_final.json
gse117606_paired_results.json enrichr_results.json consensus_all.csv consensus_core.csv
consensus_summary.json _bmc_references.json""".split()

# 三个额外队列的 DEG 结果（Stage A 输入之一）
DEG3 = ["GSE37364_DEG_Carcinoma_vs_Adenoma.csv",
        "GSE20916macro_DEG_Carcinoma_vs_Adenoma.csv",
        "GSE41657_DEG_Carcinoma_vs_Adenoma.csv"]

# 上游 R 流程产物：原绝对路径 → 仓库内相对路径
UPSTREAM = [
    ("D:/2/01.limma/DEG_all_Tumor_vs_Polyp.csv", "upstream_r/01.limma"),
    ("D:/2/02.WGCNA/WGCNA_module_genes.csv", "upstream_r/02.WGCNA"),
    ("D:/2/02.WGCNA/WGCNA_top_module_genes.csv", "upstream_r/02.WGCNA"),
    ("D:/2/03.enrich/gene_Tumor_vs_Polyp.csv", "upstream_r/03.enrich"),
    ("D:/2/03.enrich/Hub_Tumor_vs_Polyp.csv", "upstream_r/03.enrich"),
]

TREES = ["figs", "suppl"]
DOCS = ["REPRODUCIBILITY.md"]

# 树状目录里**不复制**的文件。
# 补充表索引 HTML 是**内部**文档：正文为中文，且通篇用内部章节编号（R2 ¶3、Methods M12），
# 放进公开仓库既不是英文、又泄露内部脚手架编号。它在工作区里保留，只是不进仓库。
TREE_EXCLUDE = {"Supplementary_Tables_S1-S5索引.html"}

# 仓库专属文件（README / LICENSE / 依赖清单 / gitignore），源放在 _repo_assets/ 下，
# 这样本脚本可以整体重跑而不会丢掉它们。
ASSETS_DIR = "_repo_assets"
ASSETS = ["README.md", "LICENSE", "REQUIREMENTS.txt", ".gitignore",
          ".zenodo.json", "RELEASE_NOTES_v1.0.0.md"]

copied, missing, skipped = [], [], []


def put(rel_src, rel_dst):
    s = os.path.join(SRC, rel_src)
    d = os.path.join(DST, rel_dst)
    if not os.path.exists(s):
        missing.append(rel_src)
        return
    os.makedirs(os.path.dirname(d), exist_ok=True)
    shutil.copy2(s, d)
    copied.append(rel_dst)


# ---------------------------------------------------------------- 清空重建
# 只清内容、保留 `.git/`：否则每次重跑都会把本地提交历史（以及已配置的 remote）删掉。
if os.path.isdir(DST):
    assert os.path.basename(DST) == "repo", "拒绝删除非 repo 目录"
    for entry in os.listdir(DST):
        if entry == ".git":
            continue
        p = os.path.join(DST, entry)
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
else:
    os.makedirs(DST)

for name in SCRIPTS:
    put(name + ".py", name + ".py")
for name in MIDDLE:
    put(name, name)
for name in DEG3:
    put("02.result/" + name, "02.result/" + name)
for abs_src, rel_dir in UPSTREAM:
    if not os.path.exists(abs_src):
        missing.append(abs_src)
        continue
    target = os.path.join(DST, rel_dir, os.path.basename(abs_src))
    os.makedirs(os.path.dirname(target), exist_ok=True)
    shutil.copy2(abs_src, target)
    copied.append(f"{rel_dir}/{os.path.basename(abs_src)}")
for doc in DOCS:
    put(doc, doc)
for a in ASSETS:
    put(os.path.join(ASSETS_DIR, a), a)


# ---------------------------------------------------------------- 占位符填充
# ★ 必须在**打包阶段**做，不能只交给 push_github：否则本脚本每次重跑都会拿
#   _repo_assets/ 里的原版把 push_github 已填好的 URL 打回占位符（实测发生过，
#   结果公开 README 的 "Source repository" 一格长期显示占位符）。
#   占位符在 markdown 里写作反斜杠转义形式（\<REPO_URL\>），故两种写法都替换。
#
# ★ 必须覆盖**所有**带占位符的文本资产，不能只管 README：
#   RELEASE_NOTES_v1.0.0.md 里有同一句 `git clone <REPO_URL>`，而它正是贴进
#   GitHub release 正文的那段文字 —— 只填 README 会让 release 页面公开发布
#   一个裸占位符（实测漏网，由 preflight_release.py 兜住）。
FILL_FILES = ["README.md", "RELEASE_NOTES_v1.0.0.md"]


def _substitutions():
    vals = {}
    ca = os.path.join(SRC, "code_availability.json")
    if os.path.exists(ca):
        try:
            vals = json.load(open(ca, encoding="utf-8"))
        except Exception:
            vals = {}
    doi = vals.get("doi")
    return {
        "REPO_URL": vals.get("repo_url") or "not pushed yet",
        "ARCHIVED_DOI": ("https://doi.org/" + doi) if doi
                        else "not yet minted (Zenodo assigns one to each tagged release)",
    }


def _fill_placeholders():
    sub = _substitutions()
    filled = []
    for rel in FILL_FILES:
        path = os.path.join(DST, rel)
        if not os.path.exists(path):
            continue
        s = open(path, encoding="utf-8").read()
        orig = s
        for key, val in sub.items():
            s = s.replace("\\<%s\\>" % key, val).replace("<%s>" % key, val)
        if s != orig:
            # ★ newline="" 是必须的：Windows 文本模式写回会把 LF 全变成 CRLF。
            #   真源 _repo_assets/ 里全是 LF，一旦被转成 CRLF，改一行就会在 git diff
            #   里显示成整文件重写（实测 RELEASE_NOTES 一行 URL 修复 = 116 行 diff），
            #   把真正的改动淹掉。
            open(path, "w", encoding="utf-8", newline="").write(s)
            filled.append(rel)
    return filled


_filled = _fill_placeholders()

n_tree = {}
for t in TREES:
    src_t = os.path.join(SRC, t)
    if not os.path.isdir(src_t):
        missing.append(t + "/")
        continue
    dst_t = os.path.join(DST, t)
    os.makedirs(dst_t, exist_ok=True)
    n = 0
    for f in sorted(os.listdir(src_t)):
        sp = os.path.join(src_t, f)
        if os.path.isfile(sp):
            if f in TREE_EXCLUDE:
                skipped.append(f"{t}/{f}")
                continue
            shutil.copy2(sp, os.path.join(dst_t, f))
            copied.append(f"{t}/{f}")
            n += 1
    n_tree[t] = n

# ---------------------------------------------------------------- 报告
def _size(root):
    """统计复制内容的大小，**排除 .git** —— 否则重跑时数字会混入版本库体积。"""
    return sum(os.path.getsize(os.path.join(r, f))
               for r, _ds, fs in os.walk(root)
               if ".git" not in r.split(os.sep) for f in fs)


tot = _size(DST)
print(f"仓库目录 : {DST}")
print(f"复制文件 : {len(copied)} 个 | 合计 {tot / 1024 / 1024:.2f} MB")
print(f"  · 脚本 {len(SCRIPTS)} 个")
print(f"  · 中间产物 {len(MIDDLE)} 个")
print(f"  · 三队列 DEG {len(DEG3)} 个")
print(f"  · 上游 R 产物 {len(UPSTREAM)} 个 -> upstream_r/")
for t, n in n_tree.items():
    print(f"  · {t}/ {n} 个")
if skipped:
    print(f"  · 主动排除 {len(skipped)} 个（内部文档，不进公开仓库）：")
    for s in skipped:
        print(f"      {s}")
if missing:
    print("\n[缺失]")
    for m in missing:
        print("  MISS", m)

with open(os.path.join(DST, "_MANIFEST.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(sorted(copied)) + "\n")
print("\n已写出 _MANIFEST.txt")
if _filled:
    print("已填充占位符：" + "、".join(_filled))
else:
    print("注意：没有文件被填充占位符（检查 _repo_assets/ 里是否还有 <REPO_URL>）")
