# -*- coding: utf-8 -*-
"""Figure S3：GSE303070 空间图谱的取样偏差与配对组成读数（R5 用）。
★ 英文化：`order` 里的 niche 名是**数据键**（用于 fovT.get(lab)），保持原样；
  只有 `labels` 的**值**是上屏文字，直接改英文。保存走 fig_en.finish() 守卫。
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fig_en

fig_en.use_english_font()

FIG = "figs"
os.makedirs(FIG, exist_ok=True)

NZ = json.load(open("sc303070_niche.json", encoding="utf-8"))
CP = json.load(open("sc303070_composition_final.json", encoding="utf-8"))

fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))

# ---------- A: FOV niche 构成（T vs N）----------
ax = axes[0]
fovT = NZ["fov_by_niche"]["T"]; fovN = NZ["fov_by_niche"]["N"]
order = ["Epithelial mass", "SMC epithelial rich", "Normal like", "Normal like 2",
         "Stroma", "Vasculature", "Macrophage rich", "Granulocyte rich", "TLS", "QC_fail"]
labels = {"Epithelial mass": "epithelial nest", "SMC epithelial rich": "epithelium-enriched/SMC",
          "Normal like": "normal-like 1", "Normal like 2": "normal-like 2", "Stroma": "stroma",
          "Vasculature": "vasculature", "Macrophage rich": "macrophage-rich",
          "Granulocyte rich": "granulocyte-rich", "TLS": "lymphoid (TLS)",
          "QC_fail": "QC-failed"}
colors = ["#1d4ed8", "#3b82f6", "#0d9488", "#14b8a6", "#a16207", "#ca8a04",
          "#be123c", "#e11d48", "#7c3aed", "#64748b"]
tT = sum(fovT.values()); tN = sum(fovN.values())
y = [0, 1]
left_b = np.zeros(2)
for lab, c in zip(order, colors):
    v = np.array([100 * fovT.get(lab, 0) / tT, 100 * fovN.get(lab, 0) / tN])
    ax.barh(y, v, left=left_b, color=c, edgecolor="white", height=0.40, label=labels[lab])
    for i, (vv, lb) in enumerate(zip(v, left_b)):
        if vv >= 4.5:
            ax.text(lb + vv / 2, y[i], "%.0f%%" % vv, ha="center", va="center",
                    color="white", fontsize=8.5, fontweight="bold")
    left_b += v
ax.set_yticks(y); ax.set_yticklabels(["tumour\n(%d FOV)" % tT, "normal\n(%d FOV)" % tN], fontsize=11)
ax.set_xlabel("FOV composition (%)", fontsize=11)
ax.set_title("A  Field-of-view sampling differs systematically between tumour and normal samples\n(28% of normal FOVs are QC-failed and 16% are TLS)",
             fontsize=11.5, fontweight="bold")
# ★ 图例不能放 "lower right"：会盖住 normal 那根条上的分段百分比标签。
#   改放到两根条之间的空白带（y=0 与 y=1 之间），并把条压薄留出空间。
ax.legend(fontsize=7.5, ncol=4, loc="center", framealpha=.95)
ax.set_xlim(0, 100); ax.invert_yaxis()

# ---------- B: 配对供体上皮占比 ----------
ax = axes[1]
ps = CP["per_sample"]
paired = {}
for k, v in ps.items():
    d, t = k.split("|")[:2]
    paired.setdefault(d, {})[t] = v
dons = sorted([d for d, v in paired.items() if "T" in v and "N" in v])
xs = np.arange(len(dons)); w = 0.36
et = [paired[d]["T"]["epi"] * 100 for d in dons]
en = [paired[d]["N"]["epi"] * 100 for d in dons]
ax.bar(xs - w / 2, en, w, color="#0d9488", label="normal (matched)", alpha=.9)
ax.bar(xs + w / 2, et, w, color="#be123c", label="tumour", alpha=.9)
for i, d in enumerate(dons):
    ax.plot([xs[i] - w / 2, xs[i] + w / 2], [en[i], et[i]], color="#334155", lw=1.1, ls="--", zorder=5)
    dv = et[i] - en[i]
    ax.text(xs[i], max(et[i], en[i]) + 2.2, "%+.1f" % dv, ha="center", fontsize=9,
            color="#be123c" if dv > 0 else "#0d9488", fontweight="bold")
ax.set_xticks(xs)
ax.set_xticklabels([d.replace("Donor", "D") for d in dons], fontsize=10)
ax.set_ylabel("Epithelial cell fraction (%) (QC = High)", fontsize=11)
ax.set_ylim(0, 82)
ax.axhline(0, color="#cbd5e1")
md = CP["paired"]["median_d_epi"] * 100
ax.set_title(fig_en.wrap(
    "B  Matched donors: the tumour epithelial fraction is higher (median Δ %+.1f, concordant 4/5)\n"
    "→ opposite to the naive expectation of tumour de-epithelialisation, indicating a sampling effect"
    % md, 6.6, 11.5),
    fontsize=11.5, fontweight="bold")
ax.legend(fontsize=9)

fig.suptitle(fig_en.wrap(
    "Figure S3  GSE303070 (CosMx spatial atlas) cannot serve as an unbiased composition "
    "baseline: FOV sampling is strongly correlated with sample type", 14.0, 12.5),
    fontsize=12.5, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = os.path.join(FIG, "figS3_GSE303070_sampling_bias.png")
fig_en.finish(fig, out, dpi=300)
