# -*- coding: utf-8 -*-
"""图 S3：GSE303070 空间图谱的取样偏差与配对组成读数（R5 用）。"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

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
labels = {"Epithelial mass": "上皮团", "SMC epithelial rich": "上皮富集/平滑肌",
          "Normal like": "正常样 1", "Normal like 2": "正常样 2", "Stroma": "间质",
          "Vasculature": "血管", "Macrophage rich": "巨噬富集",
          "Granulocyte rich": "粒细胞富集", "TLS": "淋巴聚集(TLS)", "QC_fail": "QC 不合格"}
colors = ["#1d4ed8", "#3b82f6", "#0d9488", "#14b8a6", "#a16207", "#ca8a04",
          "#be123c", "#e11d48", "#7c3aed", "#64748b"]
tT = sum(fovT.values()); tN = sum(fovN.values())
y = [0, 1]
left_b = np.zeros(2)
for lab, c in zip(order, colors):
    v = np.array([100 * fovT.get(lab, 0) / tT, 100 * fovN.get(lab, 0) / tN])
    ax.barh(y, v, left=left_b, color=c, edgecolor="white", height=0.55, label=labels[lab])
    for i, (vv, lb) in enumerate(zip(v, left_b)):
        if vv >= 4.5:
            ax.text(lb + vv / 2, y[i], "%.0f%%" % vv, ha="center", va="center",
                    color="white", fontsize=8.5, fontweight="bold")
    left_b += v
ax.set_yticks(y); ax.set_yticklabels(["肿瘤\n(%d FOV)" % tT, "正常\n(%d FOV)" % tN], fontsize=11)
ax.set_xlabel("视野（FOV）构成占比 %", fontsize=11)
ax.set_title("A  视野取材类型在肿瘤/正常样本间系统性不同\n（正常样本中 28% 为 QC 不合格、16% 为 TLS）",
             fontsize=11.5, fontweight="bold")
ax.legend(fontsize=8, ncol=2, loc="lower right", framealpha=.95)
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
ax.bar(xs - w / 2, en, w, color="#0d9488", label="正常（同源）", alpha=.9)
ax.bar(xs + w / 2, et, w, color="#be123c", label="肿瘤", alpha=.9)
for i, d in enumerate(dons):
    ax.plot([xs[i] - w / 2, xs[i] + w / 2], [en[i], et[i]], color="#334155", lw=1.1, ls="--", zorder=5)
    dv = et[i] - en[i]
    ax.text(xs[i], max(et[i], en[i]) + 2.2, "%+.1f" % dv, ha="center", fontsize=9,
            color="#be123c" if dv > 0 else "#0d9488", fontweight="bold")
ax.set_xticks(xs)
ax.set_xticklabels([d.replace("Donor", "D") for d in dons], fontsize=10)
ax.set_ylabel("上皮细胞占比 %（QC=High）", fontsize=11)
ax.set_ylim(0, 82)
ax.axhline(0, color="#cbd5e1")
md = CP["paired"]["median_d_epi"] * 100
ax.set_title("B  配对供体：肿瘤上皮占比反而更高（中位 Δ %+.1f，同向 4/5）\n"
             "→ 与「肿瘤去上皮化」的朴素预期相反，提示受取材影响" % md,
             fontsize=11.5, fontweight="bold")
ax.legend(fontsize=9)

fig.suptitle("图 S3  GSE303070（CosMx 空间图谱）无法充当无偏的组成基准：视野取材与样本类型强相关",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(FIG, "figS3_GSE303070取样偏差.png"), dpi=300)
print("saved", os.path.join(FIG, "figS3_GSE303070取样偏差.png"))
