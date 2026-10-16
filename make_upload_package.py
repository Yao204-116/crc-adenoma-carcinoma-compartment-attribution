# -*- coding: utf-8 -*-
"""生成 BMC 上传包 upload/。

只做"复制 + 规范命名"，不改动、不删除任何原文件（figs/ 与 suppl/ 原样保留）。
命名依据 BMC 官网规定：附加文件须命名 "Additional file n"；Figures 须按正文首见序编号上传。
用 anaconda 的 python 运行（需要 pandas + openpyxl）。
"""
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
    ("figs/fig1_study_design.png",        "Figure_1.png"),
    ("figs/fig2_56组配对同向比例.png",      "Figure_2.png"),
    ("figs/fig3_2x2分解_全组织vsLCM.png",  "Figure_3.png"),
    ("figs/fig4_逐基因归属散点.png",       "Figure_4.png"),
    ("figs/fig5_面板队列AUC热图.png",      "Figure_5.png"),
    ("figs/fig6_三层模式_腺瘤峰.png",       "Figure_6.png"),
    ("figs/fig7_单细胞P×U分解.png",        "Figure_7.png"),
    ("figs/fig8_DEFA腺瘤峰dot.png",        "Figure_8.png"),
]

# 附加文件编号 → 上传文件名 + 源（表：[(csv, 工作表名)]；图：png 路径）
# 编号与 build_bmc_manuscript.py 按正文首见序自动算出的号一致（verify_bmc.py 已核 1..8）。
AF = {
    1: ("Additional_file_1.xlsx", [("suppl/S1_313gene_attribution.csv", "Table S1")]),
    2: ("Additional_file_2.xlsx", [("suppl/S2_panel_gene_membership.csv", "Table S2")]),
    3: ("Additional_file_3.xlsx", [("suppl/S4a_sc161277_sample_QC.csv", "Table S4a"),
                                   ("suppl/S4b_sc161277_decomposition.csv", "Table S4b")]),
    4: ("Additional_file_4.xlsx", [("suppl/S5_four_cohort_performance.csv", "Table S5")]),
    5: ("Additional_file_5.png",  "figs/figS1_漂移诊断.png"),
    6: ("Additional_file_6.png",  "figs/figS2_双口径筛选矩阵.png"),
    7: ("Additional_file_7.xlsx", [("suppl/S3_panel_pattern_matrix.csv", "Table S3a"),
                                   ("suppl/S3b_paneth_sensitivity.csv", "Table S3b")]),
    8: ("Additional_file_8.png",  "figs/figS3_GSE303070取样偏差.png"),
}

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
for n, (dst, src) in sorted(AF.items()):
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

# ---- 上传清单 ----
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
          "- Main figure numbers follow **first-mention order** in the manuscript, so the",
          "  source files `fig2_*` and `fig6_*` are swapped here (Fig. 2 ↔ Fig. 6).",
          "- Additional file numbers also follow first-mention order, so they do **not** match",
          "  the internal `Sn` labels (e.g. Additional file 3 = Table S4, Additional file 7 = Table S3).",
          "- Figures are 300 dpi PNG; BMC accepts PNG at ≤10 MB per file — no conversion needed.",
          "- Tables are delivered as .xlsx (BMC accepts .xls / .csv for tabular additional files)."]
open(os.path.join(UP, "UPLOAD_MANIFEST.md"), "w", encoding="utf-8").write("\n".join(lines))

print(f"upload/ 生成 {len(done)} 个文件")
for d, s, n in done:
    print(f"  {d:28} {s:>12}  {n}")
print("missing:", missing if missing else "无")
