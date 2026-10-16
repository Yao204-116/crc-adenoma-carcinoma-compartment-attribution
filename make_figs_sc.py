# -*- coding: utf-8 -*-
"""图 7（P×U 分解）与 图 8（DEFA5/6 dot plot）——GSE161277 单细胞章 v2。

v2 变更：
  - 适配 sc161277_pairs.json 的合并行（defa，P3 双腺瘤加权合并）与逐样本行（defa_persample）
  - 图 7A 加入 Patient3（低解析度，斜纹柱 + 脚注）
  - 图 8 用逐样本行，P3 两个腺瘤样本各自成行（不合并）
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)
for cand in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC"]:
    if any(cand.lower() in f.name.lower() for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = cand
        break

P = json.load(open(os.path.join(HERE, "sc161277_pairs.json"), encoding="utf-8"))
labels = {"normal": "正常", "adenoma": "腺瘤", "carcinoma": "癌", "para-cancer": "癌旁", "blood": "血"}

# ---------- 图 7 ----------
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
           edgecolor="white", label="组成项 Δlog P" if i == 0 else None)
    ax.bar([i + w / 2], [d["dlogU"]], w, color="#0d9488", alpha=.85, hatch=hatch,
           edgecolor="white", label="状态项 Δlog U" if i == 0 else None)
    ax.plot([i], [d["dlogT"]], marker="D", color="#e11d48", ms=8, zorder=5,
            label="Δlog T（LCM 视角）" if i == 0 else None)
    ax.text(i - w / 2, d["dlogP"] + (0.06 if d["dlogP"] >= 0 else -0.14), f"{d['dlogP']:+.2f}",
            ha="center", fontsize=10)
    ax.text(i + w / 2, d["dlogU"] + (0.06 if d["dlogU"] >= 0 else -0.14), f"{d['dlogU']:+.2f}",
            ha="center", fontsize=10, fontweight="bold")
ax.axhline(0, color="#334155", lw=1)
ax.set_xticks(range(len(dec)))
ax.set_xticklabels(["%s\n(上皮占比 %.1f%%→%.1f%%)" % (d["patient"], d["P_a"] * 100, d["P_b"] * 100)
                    for d in dec], fontsize=9.5)
ax.set_xlabel("斜纹柱 = 上皮占比 <4%，分解解析度不足（仅供方向参考）", fontsize=8.5, color="#b45309")
ax.set_ylabel("潘氏签名 Δlog（正常→腺瘤）")
ax.set_title("A  正常→腺瘤：潘氏峰 = 组成 + 状态", fontsize=12, fontweight="bold")
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
ax.set_ylabel("DEFA6 上皮内表达（U 层）")
ax.set_title("B  DEFA6：腺瘤特异升高、癌期关闭（3/3 患者）", fontsize=12, fontweight="bold")
ax.legend(fontsize=9)
ax.margins(y=0.18)      # 留白避免柱顶数值标签被裁切
fig.suptitle("GSE161277 单细胞分解：腺瘤期潘氏峰的结构（三患者配对 normal→腺瘤→癌）",
             fontsize=13.5, fontweight="bold", y=1.00)
fig.text(0.995, 0.01, "T=P×U：Δlog T = Δlog P（组成）+ Δlog U（状态）；DEFA6 为纯潘氏基因、无免疫污染",
         ha="right", fontsize=8.5, color="#64748b")
fig.tight_layout(rect=(0, 0.03, 1, 0.97))
fig.savefig(os.path.join(FIG, "fig7_单细胞P×U分解.png"), dpi=300)
print("fig7 saved")

# ---------- 图 8：逐样本 DEFA5/DEFA6 dot（U 层，P3 两腺瘤各自成行） ----------
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

fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), sharey=False)
for k, (gene, ttl) in enumerate([("DEFA5_U", "DEFA5"), ("DEFA6_U", "DEFA6")]):
    ax = axes[k]
    for i, r in enumerate(order):
        v = r[gene]
        col = "#e11d48" if r["stage"] == "adenoma" else ("#94a3b8" if r["stage"] == "blood" else "#2563eb")
        ax.scatter(v, i, s=max(28, r["n_epi"] / 14), c=col, alpha=.85,
                   edgecolors="white", linewidths=1.2, zorder=3)
        ax.text(v, i - 0.36, f"{v:.3f}", fontsize=7.5, ha="center", color="#475569")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(ylab, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color="#cbd5e1", lw=1)
    ax.set_xlabel(f"{gene} 上皮内表达（U 层）")
    ax.set_title(f"{ttl}：腺瘤（红）vs 正常/癌（蓝）", fontsize=11, fontweight="bold")
from matplotlib.lines import Line2D
leg = [Line2D([], [], marker="o", ls="", color="#e11d48", label="腺瘤"),
       Line2D([], [], marker="o", ls="", color="#2563eb", label="正常/癌/癌旁"),
       Line2D([], [], marker="o", ls="", color="#94a3b8", label="血（阴性）")]
axes[0].legend(handles=leg, fontsize=8.5, loc="lower right")
fig.suptitle("纯潘氏基因在上皮细胞内的表达（GSE161277，逐样本；点面积 ∝ 上皮细胞数）",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(FIG, "fig8_DEFA腺瘤峰dot.png"), dpi=300)
print("fig8 saved")
