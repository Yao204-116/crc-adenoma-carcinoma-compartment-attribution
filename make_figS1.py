# -*- coding: utf-8 -*-
"""Figure S1：全局漂移诊断（R2 / Methods M3 的支撑图）。
A  三队列全基因组 Δ(癌−腺瘤) 与上调基因占比——极不均衡
B  GSE15960 稳定基因探针诊断——偏移为技术性，校正只吸收部分落差
数据：manuscript_facts.json global_drift ／ drift_diag.json
★ 英文化：数据键（队列代码 / stage 代码）不动，只翻上屏文字；保存走 fig_en.finish() 守卫。
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

GD = json.load(open("manuscript_facts.json", encoding="utf-8"))["global_drift"]
DD = json.load(open("drift_diag.json", encoding="utf-8"))["GSE15960"]

fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6),
                         gridspec_kw={"width_ratios": [1.2, 1]})

# ---------- A：三队列全基因组漂移 ----------
ax = axes[0]
COH = ["GSE4183", "GSE117606", "GSE15960"]
LAB = {"GSE4183": "GSE4183\nwhole tissue·Affymetrix",
       "GSE117606": "GSE117606\nwhole tissue·Agilent",
       "GSE15960": "GSE15960\nLCM pure epithelium·Affymetrix"}
COL = ["#1d4ed8", "#3b82f6", "#be123c"]
d = [GD[c]["global_delta"] for c in COH]
sd = [GD[c]["global_sd"] for c in COH]
fu = [GD[c]["frac_up"] for c in COH]
x = np.arange(3)
ax.bar(x, d, 0.46, color=COL, alpha=.92, edgecolor="white",
       label="genome-wide Δ (carcinoma−adenoma, z units)\namber ↑ = up-regulated gene fraction (%)")
for i in range(3):
    y0 = min(d[i], 0)
    # ★ 英文标签比中文长得多：首/末组整组外移 0.30 并改左/右对齐，否则相邻两组会互相压叠
    #   （守卫的压叠检测实测：'…42.9%' × '…37.4%' 重叠 35%）。xlim 约 (−0.35, 2.35)，0.30 是安全上限。
    xa = x[i] - 0.30 if i == 0 else (x[i] + 0.30 if i == 2 else x[i])
    ha = "left" if i == 0 else ("right" if i == 2 else "center")
    ax.text(xa, y0 - 0.028, "Δ %.3f | SD %.3f" % (d[i], sd[i]),
            ha=ha, va="top", fontsize=10, fontweight="bold", color="#334155")
    ax.text(xa, y0 - 0.115, "↑ %.1f%%" % fu[i],
            ha=ha, va="top", fontsize=10, fontweight="bold", color="#a16207")
ax.axhline(0, color="#94a3b8", lw=1)
ax.set_xticks(x); ax.set_xticklabels([LAB[c] for c in COH], fontsize=10)
ax.set_ylabel("Genome-wide inter-group offset Δ (z units)", fontsize=11)
ax.set_ylim(-0.63, 0.15)
ax.set_yticks([0, -0.1, -0.2, -0.3, -0.4, -0.5])
ax.legend(fontsize=9, loc="lower left", framealpha=.95)
ax.set_title(fig_en.wrap(
    "A  Group offsets are highly imbalanced across the three cohorts: LCM Δ = −0.400 "
    "(others ≈ −0.06)\n"
    "→ rule out a global intensity offset between sample groups before invoking biology",
    7.0, 11.5), fontsize=11.5, fontweight="bold")

# ---------- B：稳定基因探针诊断 ----------
ax = axes[1]
G = ["normal", "adenoma", "carcinoma"]
GN = ["Normal", "Adenoma", "Carcinoma"]
raw = [DD["group_glev"][g] for g in G]
adj = [DD["stable_glev_mean"][g] for g in G]
x = np.arange(3); w = 0.36
ax.bar(x - w / 2, raw, w, color="#64748b", alpha=.9, label="uncorrected (group means)")
ax.bar(x + w / 2, adj, w, color="#0d9488", alpha=.9, label="residual after sample-level correction")
for i in range(3):
    ax.text(x[i] - w / 2, raw[i] + (0.012 if raw[i] >= 0 else -0.012),
            "%+.3f" % raw[i], ha="center",
            va="bottom" if raw[i] >= 0 else "top", fontsize=9.5)
    ax.text(x[i] + w / 2, adj[i] + (0.012 if adj[i] >= 0 else -0.012),
            "%+.3f" % adj[i], ha="center",
            va="bottom" if adj[i] >= 0 else "top", fontsize=9.5, color="#0f766e")
ax.axhline(0, color="#94a3b8", lw=1)
ax.text(1.0, -0.30, "uncorrected adenoma−carcinoma Δ = −0.400\nresidual after correction Δ = +0.171",
        ha="center", va="center", fontsize=10, fontweight="bold", color="#334155")
ax.set_xticks(x); ax.set_xticklabels(GN, fontsize=11)
ax.set_ylabel("Group means over stable gene probes (z units)", fontsize=11)
ax.set_ylim(-0.35, 0.34)
ax.legend(fontsize=9, loc="upper right", framealpha=.95)
ax.set_title(fig_en.wrap(
    "B  Stable gene probes (high expression, smallest between-sample SD) reproduce −0.400 "
    "on the raw scale\n"
    "→ the offset is technical; the residual is +0.171, hence both calibers are reported throughout",
    5.9, 11.5), fontsize=11.5, fontweight="bold")

fig.suptitle(fig_en.wrap(
    "Figure S1  Global drift diagnostics: the −0.400 group offset in the LCM cohort is "
    "technical and is not fully corrected", 14.0, 13.0),
    fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
out = os.path.join(FIG, "figS1_drift_diagnostics.png")
fig_en.finish(fig, out, dpi=300)
print("A  global drift:", {c: (round(GD[c]["global_delta"], 3), round(GD[c]["frac_up"], 1)) for c in COH})
print("B  uncorrected:", {g: round(DD["group_glev"][g], 3) for g in G},
      " corrected:", {g: round(DD["stable_glev_mean"][g], 3) for g in G})
