# -*- coding: utf-8 -*-
"""Figure 7（P×U 分解）与 Figure 8（DEFA5/6 dot plot）——GSE161277 单细胞章 v2。

v2 变更：
  - 适配 sc161277_pairs.json 的合并行（defa，P3 双腺瘤加权合并）与逐样本行（defa_persample）
  - Figure 7A 加入 Patient3（低解析度，斜纹柱 + 脚注）
  - Figure 8 用逐样本行，P3 两个腺瘤样本各自成行（不合并）

★ 英文化：`labels` 字典的**值**只用于上屏（键仍是数据里的 stage 代码），所以直接改英文是安全的；
  两处 f-string（xlabel / title）求值后才含中文，改成显示名 `ttl` 拼接。
  保存统一走 `fig_en.finish()` 守卫。
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fig_en

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)
fig_en.use_english_font()

P = json.load(open(os.path.join(HERE, "sc161277_pairs.json"), encoding="utf-8"))
# 键 = 数据里的 stage 代码；值 = 上屏文字（英文）
labels = {"normal": "Normal", "adenoma": "Adenoma", "carcinoma": "Carcinoma",
          "para-cancer": "para-cancer", "blood": "Blood"}

# ---------- Figure 7 ----------
dec = [d for d in P["decomp"] if d["sig"] == "paneth" and d["a"] == "normal"]
dec = sorted(dec, key=lambda d: d["patient"])
patients = [d["patient"] for d in dec]
LOWRES = {"Patient3"}          # 上皮占比 <4%，分解解析度不足
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

ax = axes[0]
w = 0.36
for i, d in enumerate(dec):
    low = d["patient"] in LOWRES
    hatch = "///" if low else None
    ax.bar([i - w / 2], [d["dlogP"]], w, color="#7c3aed", alpha=.85, hatch=hatch,
           edgecolor="white", label="composition term Δlog P" if i == 0 else None)
    ax.bar([i + w / 2], [d["dlogU"]], w, color="#0d9488", alpha=.85, hatch=hatch,
           edgecolor="white", label="state term Δlog U" if i == 0 else None)
    ax.plot([i], [d["dlogT"]], marker="D", color="#e11d48", ms=8, zorder=5,
            label="Δlog T (LCM view)" if i == 0 else None)
    ax.text(i - w / 2, d["dlogP"] + (0.06 if d["dlogP"] >= 0 else -0.14), f"{d['dlogP']:+.2f}",
            ha="center", fontsize=10)
    ax.text(i + w / 2, d["dlogU"] + (0.06 if d["dlogU"] >= 0 else -0.14), f"{d['dlogU']:+.2f}",
            ha="center", fontsize=10, fontweight="bold")
ax.axhline(0, color="#334155", lw=1)
ax.set_xticks(range(len(dec)))
# ★ 不能再用「PatientN\n(epithelial fraction a%→b%)」：英文串太长，三组标签会互相压叠（目视复查发现）。
#   拆成三行短串，并在 xlabel 里定义 "epithelium"。
ax.set_xticklabels(["%s\nepithelium\n%.1f%%→%.1f%%" % (d["patient"], d["P_a"] * 100, d["P_b"] * 100)
                    for d in dec], fontsize=9)
ax.set_xlabel("hatched bars = epithelial-cell fraction <4%, decomposition under-resolved (direction only)\n"
              "epithelium = epithelial-cell fraction (%) of the sample",
              fontsize=8.5, color="#b45309")
ax.set_ylabel("Paneth signature Δlog (normal→adenoma)")
ax.set_title(fig_en.wrap("A  normal→adenoma: the Paneth peak = composition + state", 5.5, 12),
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper right")

# 右：DEFA6 分期（合并行：P3 两个腺瘤加权合并）
ax = axes[1]
M = {(r["patient"], r["stage"]): r for r in P["defa"]}
stages = ["normal", "adenoma", "carcinoma"]
pats = ["Patient1", "Patient2", "Patient3"]
colors = {"Patient1": "#2563eb", "Patient2": "#0d9488", "Patient3": "#d97706"}
w = 0.26
for j, p in enumerate(pats):
    vals = [M[(p, s)]["DEFA6_U"] if (p, s) in M else 0 for s in stages]
    xs = [x + (j - 1) * w for x in range(len(stages))]
    ax.bar(xs, vals, w, color=colors[p], alpha=.85, label=p)
    for x, v in zip(xs, vals):
        ax.text(x, v + 0.007, f"{v:.3f}", ha="center", fontsize=8.5, rotation=90)
ax.set_xticks(range(len(stages)))
ax.set_xticklabels([labels[s] for s in stages], fontsize=11)
ax.set_ylabel("DEFA6 intra-epithelial expression (U term)")
ax.set_title(fig_en.wrap("B  DEFA6: adenoma-specific rise, switched off at carcinoma (3/3 patients)",
                         5.5, 11.5),
             fontsize=11.5, fontweight="bold")
ax.legend(fontsize=9)
ax.margins(y=0.18)      # 留白避免柱顶数值标签被裁切
fig.suptitle(fig_en.wrap(
    "GSE161277 single-cell decomposition: structure of the adenoma-stage Paneth peak "
    "(three patients, paired normal→adenoma→carcinoma)", 12.4, 13.0),
    fontsize=13.0, fontweight="bold", y=0.99)
fig.text(0.995, 0.01, "T = P × U: Δlog T = Δlog P (composition) + Δlog U (state); DEFA6 is a pure Paneth gene with no immune contamination",
         ha="right", fontsize=8.5, color="#64748b")
fig.tight_layout(rect=(0, 0.03, 1, 0.94))
fig_en.finish(fig, os.path.join(FIG, "fig7_single_cell_PxU_decomposition.png"), dpi=300)

# ---------- Figure 8：逐样本 DEFA5/DEFA6 dot（U 层，P3 两腺瘤各自成行） ----------
PS = P["defa_persample"]
key = {(r["patient"], r["stage"], r["gsm"]): r for r in PS}
order = []
for p in ["Patient1", "Patient2", "Patient3"]:
    for s in ["normal", "para-cancer", "adenoma", "carcinoma"]:
        for r in PS:
            if r["patient"] == p and r["stage"] == s:
                order.append(r)
for r in PS:
    if r["patient"] == "Patient0":
        order.append(r)
for r in PS:
    if r["stage"] == "blood":
        order.append(r)

dup = {}
for r in order:                       # P3 两个腺瘤区分标注
    k = (r["patient"], r["stage"])
    dup[k] = dup.get(k, 0) + 1
seen = {}
ylab = []
for r in order:
    k = (r["patient"], r["stage"])
    seen[k] = seen.get(k, 0) + 1
    suf = "" if dup[k] == 1 else (" #%d" % seen[k])
    ylab.append("%s\n%s%s" % (r["patient"], labels.get(r["stage"], r["stage"]), suf))

fig, axes = plt.subplots(1, 2, figsize=(13.5, 7.2), sharey=False)
for k, (gene, ttl) in enumerate([("DEFA5_U", "DEFA5"), ("DEFA6_U", "DEFA6")]):
    ax = axes[k]
    for i, r in enumerate(order):
        v = r[gene]
        col = "#e11d48" if r["stage"] == "adenoma" else ("#94a3b8" if r["stage"] == "blood" else "#2563eb")
        ax.scatter(v, i, s=max(28, r["n_epi"] / 14), c=col, alpha=.85,
                   edgecolors="white", linewidths=1.2, zorder=3)
        ax.text(v, i - 0.36, f"{v:.3f}", fontsize=7.5, ha="center", color="#475569")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(ylab, fontsize=7.5)      # ★ 行距紧：17 行 × 2 行标签，8 号会几乎贴在一起
    ax.invert_yaxis()
    ax.axvline(0, color="#cbd5e1", lw=1)
    ax.margins(x=0.06)          # ★ 留出右侧边距：否则靠右的数值标签会被画布裁掉
    ax.set_xlabel(f"{ttl} intra-epithelial expression (U term)")
    ax.set_title(fig_en.wrap(f"{ttl}: adenoma (red) vs normal/carcinoma (blue)", 5.8, 11),
                 fontsize=11, fontweight="bold")
from matplotlib.lines import Line2D
leg = [Line2D([], [], marker="o", ls="", color="#e11d48", label="Adenoma"),
       Line2D([], [], marker="o", ls="", color="#2563eb", label="normal/carcinoma/para-cancer"),
       Line2D([], [], marker="o", ls="", color="#94a3b8", label="Blood (negative control)")]
axes[0].legend(handles=leg, fontsize=8.5, loc="lower right")
fig.suptitle("Pure Paneth genes DEFA5 / DEFA6 within epithelial cells (GSE161277, per sample; dot area ∝ epithelial cell count)",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig_en.finish(fig, os.path.join(FIG, "fig8_DEFA6_adenoma_peak_dot.png"), dpi=300)
print("fig7 + fig8 written")
