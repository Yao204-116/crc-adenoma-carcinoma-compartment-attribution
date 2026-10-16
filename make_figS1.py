# -*- coding: utf-8 -*-
"""图 S1：全局漂移诊断（R2 / Methods M3 的支撑图）。
A  三队列全基因组 Δ(癌−腺瘤) 与上调基因占比——极不均衡
B  GSE15960 稳定基因探针诊断——偏移为技术性，校正只吸收部分落差
数据：manuscript_facts.json global_drift ／ drift_diag.json
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

FIG = "figs"
os.makedirs(FIG, exist_ok=True)

GD = json.load(open("manuscript_facts.json", encoding="utf-8"))["global_drift"]
DD = json.load(open("drift_diag.json", encoding="utf-8"))["GSE15960"]

fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6),
                         gridspec_kw={"width_ratios": [1.2, 1]})

# ---------- A：三队列全基因组漂移 ----------
ax = axes[0]
COH = ["GSE4183", "GSE117606", "GSE15960"]
LAB = {"GSE4183": "GSE4183\n全组织·Affymetrix", "GSE117606": "GSE117606\n全组织·Agilent",
       "GSE15960": "GSE15960\nLCM 纯上皮·Affymetrix"}
COL = ["#1d4ed8", "#3b82f6", "#be123c"]
d = [GD[c]["global_delta"] for c in COH]
sd = [GD[c]["global_sd"] for c in COH]
fu = [GD[c]["frac_up"] for c in COH]
x = np.arange(3)
ax.bar(x, d, 0.46, color=COL, alpha=.92, edgecolor="white",
       label="全基因组 Δ（癌−腺瘤，z 单位）")
for i in range(3):
    y0 = min(d[i], 0)
    ax.text(x[i], y0 - 0.028, "Δ %.3f   ｜   SD %.3f" % (d[i], sd[i]),
            ha="center", va="top", fontsize=10, fontweight="bold", color="#334155")
    ax.text(x[i], y0 - 0.115, "上调基因占比 %.1f%%" % fu[i],
            ha="center", va="top", fontsize=10, fontweight="bold", color="#a16207")
ax.axhline(0, color="#94a3b8", lw=1)
ax.set_xticks(x); ax.set_xticklabels([LAB[c] for c in COH], fontsize=10)
ax.set_ylabel("全基因组组间偏移 Δ（z 单位）", fontsize=11)
ax.set_ylim(-0.63, 0.15)
ax.set_yticks([0, -0.1, -0.2, -0.3, -0.4, -0.5])
ax.legend(fontsize=9, loc="lower left", framealpha=.95)
ax.set_title("A  三队列的组间偏移极不均衡：LCM 队列 Δ = −0.400（其余约 −0.06）\n"
             "→ 先排除“样本组整体强度偏移”，再谈生物学", fontsize=11.5, fontweight="bold")

# ---------- B：稳定基因探针诊断 ----------
ax = axes[1]
G = ["normal", "adenoma", "carcinoma"]
GN = ["正常", "腺瘤", "癌"]
raw = [DD["group_glev"][g] for g in G]
adj = [DD["stable_glev_mean"][g] for g in G]
x = np.arange(3); w = 0.36
ax.bar(x - w / 2, raw, w, color="#64748b", alpha=.9, label="未校正（组均值）")
ax.bar(x + w / 2, adj, w, color="#0d9488", alpha=.9, label="样本级校正后残差")
for i in range(3):
    ax.text(x[i] - w / 2, raw[i] + (0.012 if raw[i] >= 0 else -0.012),
            "%+.3f" % raw[i], ha="center",
            va="bottom" if raw[i] >= 0 else "top", fontsize=9.5)
    ax.text(x[i] + w / 2, adj[i] + (0.012 if adj[i] >= 0 else -0.012),
            "%+.3f" % adj[i], ha="center",
            va="bottom" if adj[i] >= 0 else "top", fontsize=9.5, color="#0f766e")
ax.axhline(0, color="#94a3b8", lw=1)
ax.text(1.0, -0.315, "未校正 腺瘤−癌 Δ = −0.400   ｜   校正后残差 Δ = +0.171",
        ha="center", va="center", fontsize=10.5, fontweight="bold", color="#334155")
ax.set_xticks(x); ax.set_xticklabels(GN, fontsize=11)
ax.set_ylabel("稳定基因集上的组均值（z 单位）", fontsize=11)
ax.set_ylim(-0.35, 0.34)
ax.legend(fontsize=9, loc="upper right", framealpha=.95)
ax.set_title("B  稳定基因探针（表达高、样本间 SD 最小）在原值上重现 −0.400\n"
             "→ 偏移是技术性的；校正后残差 +0.171，故全程双口径报告", fontsize=11.5, fontweight="bold")

fig.suptitle("图 S1  全局漂移诊断：LCM 队列的 −0.400 组间偏移为技术性，校正不完全",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.93))
out = os.path.join(FIG, "figS1_漂移诊断.png")
fig.savefig(out, dpi=300)
print("saved", out)
print("A: ", {c: (round(GD[c]["global_delta"], 3), round(GD[c]["frac_up"], 1)) for c in COH})
print("B 未校正:", {g: round(DD["group_glev"][g], 3) for g in G},
      " 校正:", {g: round(DD["stable_glev_mean"][g], 3) for g in G})
