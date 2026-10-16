# -*- coding: utf-8 -*-
"""把复现链打包成可直接推送到 GitHub 的仓库目录。

原则：**只复制，绝不修改原文件**（两个脚本的路径可移植化改造已在源文件完成）。
产物：`repo/`
"""
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
make_repo push_github""".split()

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

# 仓库专属文件（README / LICENSE / 依赖清单 / gitignore），源放在 _repo_assets/ 下，
# 这样本脚本可以整体重跑而不会丢掉它们。
ASSETS_DIR = "_repo_assets"
ASSETS = ["README.md", "LICENSE", "REQUIREMENTS.txt", ".gitignore",
          ".zenodo.json", "RELEASE_NOTES_v1.0.0.md"]

copied, missing = [], []


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
if missing:
    print("\n[缺失]")
    for m in missing:
        print("  MISS", m)

with open(os.path.join(DST, "_MANIFEST.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(sorted(copied)) + "\n")
print("\n已写出 _MANIFEST.txt")
