# -*- coding: utf-8 -*-
"""生成 BMC 上传包 upload/。

只做"复制 + 规范命名"，不改动、不删除任何原文件（figs/ 与 suppl/ 原样保留）。
命名依据 BMC 官网规定：附加文件须命名 "Additional file n"；Figures 须按正文首见序编号上传。
用 anaconda 的 python 运行（需要 pandas + openpyxl）。
"""
import datetime
import os
import re

import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)
UP = "upload"
os.makedirs(UP, exist_ok=True)

MS = "Manuscript_BMC_Genomics_投稿版.html"
ms = open(MS, encoding="utf-8").read()

# 主图：源文件名 → 正文图号（figs/ 里文件名已按最终图号改好，故这里是恒等映射）
MAIN_FIGS = [
    ("figs/fig1_study_design.png",                          "Figure_1.png"),
    ("figs/fig2_56_pair_concordance.png",                   "Figure_2.png"),
    ("figs/fig3_2x2_decomposition_whole_tissue_vs_LCM.png", "Figure_3.png"),
    ("figs/fig4_per_gene_attribution_scatter.png",          "Figure_4.png"),
    ("figs/fig5_panel_cohort_AUC_heatmap.png",              "Figure_5.png"),
    ("figs/fig6_three_layer_pattern_adenoma_peak.png",      "Figure_6.png"),
    ("figs/fig7_single_cell_PxU_decomposition.png",         "Figure_7.png"),
    ("figs/fig8_DEFA6_adenoma_peak_dot.png",                "Figure_8.png"),
]

# 附加文件编号 → 上传文件名 + 源（表：[(csv, 工作表名)]；图：png 路径）
# ★ 编号**绝不写死**。编号由 build_bmc_manuscript.py 按正文首见序算出，这里只是
#   「读回稿件的 n→标签」再查「标签→源文件」的登记表。
#   2026-10-20 踩过的坑：旧版把 n→文件 写死成 7=Table S3 / 8=Fig S3，
#   稿件重排后 Table S3 掉到 Figure S3 之后，结果 7 装表、8 装图，与正文所指**完全相反**。
#   这类错误构建不报错、自检也只看「1..8 齐全」，只有拿正文逐个核对才发现。
TABLE_SRC = {
    "Table S1": [("suppl/S1_313gene_attribution.csv", "Table S1")],
    "Table S2": [("suppl/S2_panel_gene_membership.csv", "Table S2")],
    "Table S3": [("suppl/S3_panel_pattern_matrix.csv", "Table S3a"),
                 ("suppl/S3b_paneth_sensitivity.csv", "Table S3b")],
    "Table S4": [("suppl/S4a_sc161277_sample_QC.csv", "Table S4a"),
                 ("suppl/S4b_sc161277_decomposition.csv", "Table S4b")],
    "Table S5": [("suppl/S5_four_cohort_performance.csv", "Table S5")],
}
FIG_SRC = {
    "Figure S1": "figs/figS1_drift_diagnostics.png",
    "Figure S2": "figs/figS2_double_caliber_screening_matrix.png",
    "Figure S3": "figs/figS3_GSE303070_sampling_bias.png",
}


def af_plan():
    """n → (上传文件名, 源)。标签全部来自稿件，源全部来自登记表，两边必须一一对应。"""
    plan = {}
    used = set()
    for n in range(1, 9):
        label = titles[n][0].strip()
        if label in TABLE_SRC:
            plan[n] = ("Additional_file_%d.xlsx" % n, TABLE_SRC[label])
        elif label in FIG_SRC:
            plan[n] = ("Additional_file_%d.png" % n, FIG_SRC[label])
        else:
            raise SystemExit("附加文件 %d 的标签 %r 不在登记表里（表=%s；图=%s）"
                             % (n, label, sorted(TABLE_SRC), sorted(FIG_SRC)))
        used.add(label)
    unused = (set(TABLE_SRC) | set(FIG_SRC)) - used
    assert not unused, "登记表里有稿件没用到的标签：%s" % sorted(unused)
    return plan


# 从稿件读回附加文件标题，写进清单（保证与正文一致）
titles = {}
for m in re.finditer(r"<li><b>Additional file (\d+)</b>: (.*?) \((.*?)\) &mdash; (.*?)</li>", ms):
    titles[int(m.group(1))] = (m.group(2), m.group(3), m.group(4))
assert sorted(titles) == list(range(1, 9)), f"稿件里只解析到 {sorted(titles)}"

missing, done = [], []

# ---- 主图 ----
for src, dst in MAIN_FIGS:
    if not os.path.exists(src):
        missing.append(src)
        continue
    with open(src, "rb") as f, open(os.path.join(UP, dst), "wb") as g:
        g.write(f.read())
    done.append((dst, f"{os.path.getsize(src):,} B", "<- " + src))

# ---- 附加文件 ----
for n, (dst, src) in sorted(af_plan().items()):
    out = os.path.join(UP, dst)
    if isinstance(src, str):                       # 图
        if not os.path.exists(src):
            missing.append(src)
            continue
        with open(src, "rb") as f, open(out, "wb") as g:
            g.write(f.read())
        note = "<- " + src
    else:                                          # 表 → 多 sheet xlsx
        sheets = {}
        for csv, sheet in src:
            if not os.path.exists(csv):
                missing.append(csv)
                continue
            sheets[sheet] = pd.read_csv(csv, encoding="utf-8-sig")
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            for sheet, df in sheets.items():
                df.to_excel(w, sheet_name=sheet, index=False)
        note = "<- " + " + ".join(c for c, _ in src)
    done.append((dst, f"{os.path.getsize(out):,} B", note))

# ---- 清掉上一轮遗留的错名文件（只移走，不删除）----
# 改了编号规则后，旧名文件会残留在 upload/ 里，造成「7 有两个版本」这种新问题。
expected = {dst for dst, _ in af_plan().values()} | {d for _, d in MAIN_FIGS}
stale_dir = os.path.join(UP, "_stale_%s" % datetime.date.today())
stale = [f for f in os.listdir(UP)
         if re.match(r"^(Additional_file_\d+\.(xlsx|png)|Figure_\d+\.png)$", f)
         and f not in expected]
if stale:
    os.makedirs(stale_dir, exist_ok=True)
    for f in stale:
        os.replace(os.path.join(UP, f), os.path.join(stale_dir, f))
    print(f"移走 {len(stale)} 个旧命名文件 → {stale_dir}/：{sorted(stale)}")

# ---- 上传清单 ----
plan = af_plan()
mis = [f"Additional file {n} = {titles[n][0]}" for n in sorted(plan)
       if not titles[n][0].endswith("S%d" % n)]
eg = "；".join(mis[:2]) if mis else "无"
lines = ["# BMC Genomics — upload package", "",
         "All files in this folder are **copies**; `figs/` and `suppl/` are untouched.", "",
         "## Main figures (upload in this order)", "",
         "| File | Size | Source |", "|---|---|---|"]
for d, s, n in done[:8]:
    lines.append(f"| `{d}` | {s} | `{n[3:]}` |")
lines += ["", "## Additional files", "",
          "| File | Title of data | Format | Size | Source |", "|---|---|---|---|---|"]
for d, s, n in done[8:]:
    num = int(re.search(r"_(\d+)\.", d).group(1))
    label, fmt, title = titles[num]
    lines.append(f"| `{d}` | {label} — {title} | {fmt} | {s} | `{n[3:]}` |")
lines += ["", "## Notes", "",
          "- Main figure numbers follow **first-mention order** in the manuscript. The source",
          "  file names under `figs/` already carry the final figure number, so the mapping",
          "  above is one-to-one (`figs/figN_*.png` → `Figure_N.png`).",
          "- Additional file numbers also follow first-mention order, so they do **not** match",
          "  the internal `Sn` labels (%s)." % eg,
          "- Figures are 300 dpi PNG; BMC accepts PNG at ≤10 MB per file — no conversion needed.",
          "- Tables are delivered as .xlsx (BMC accepts .xls / .csv for tabular additional files)."]
open(os.path.join(UP, "UPLOAD_MANIFEST.md"), "w", encoding="utf-8").write("\n".join(lines))

print(f"upload/ 生成 {len(done)} 个文件")
for d, s, n in done:
    print(f"  {d:28} {s:>12}  {n}")
print("\n附加文件编号 → 稿件里的标签（拿这个核对正文所指）：")
for n in sorted(plan):
    print("  %-20s → %s" % ("Additional file %d" % n, titles[n][0]))
print("missing:", missing if missing else "无")
