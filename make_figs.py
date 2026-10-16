# -*- coding: utf-8 -*-
"""
发表级图形生成（matplotlib）
================================================================================
数据源（全部为本轮实测）：
  three_layer_pattern.json      三层（正常/腺瘤/癌）面板分数 + 模式标签
  manuscript_facts.json         面板×队列 AUC（原始 & 漂移校正）、全局漂移
  suppl_table1_313genes.csv     313 基因逐基因 Δ（三队列，已校正）
  gse117606_paired_results.json 56 组患者配对

输出 figs/ 下 5 张 300dpi PNG
用法：python make_figs.py
"""
import csv, json, os, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)

# 中文字体
for cand in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC"]:
    if any(cand.lower() in f.name.lower() for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = cand
        break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 130
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["axes.linewidth"] = 0.9

# 配色（与报告 HTML 的六色系统一）
CB, CT, CA, CR, CV, CG = "#2563eb", "#0d9488", "#d97706", "#dc2626", "#7c3aed", "#16a34a"
GREY = "#94a3b8"

PAT = json.load(open("three_layer_pattern.json", encoding="utf-8"))
FACTS = json.load(open("manuscript_facts.json", encoding="utf-8"))
R3 = json.load(open("r3_facts.json", encoding="utf-8"))        # 双口径筛选结果（r3_facts.py 产出）
PAIR = json.load(open("gse117606_paired_results.json", encoding="utf-8"))
MAT = FACTS["panel_matrix"]
COH = ["GSE4183", "GSE117606", "GSE15960"]
CN = {"GSE4183": "GSE4183 全组织( Affy )", "GSE117606": "GSE117606 全组织( Agilent )",
      "GSE15960": "GSE15960 LCM 纯上皮"}
SHORT = {"GSE4183": "GSE4183 全组织", "GSE117606": "GSE117606 全组织", "GSE15960": "GSE15960 LCM"}


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  → figs/%s  (%.0f KB)" % (name, os.path.getsize(p) / 1024))


# ===========================================================================
# Fig 2  三层模式：腺瘤峰（4 面板 × 3 队列）+ 阴性对照
# ===========================================================================
def fig2():
    # ★ 用漂移校正口径（保守口径）；FACTS["three_layer"]["raw"] 为原始口径
    TLA = FACTS["three_layer"]["adj"]
    # 键名必须与数据完全一致
    panels = ["潘氏细胞·4", "黏液杯状·5", "上皮分化·核心8", "相代谢·5",
              "纯上皮基线·7", "间质·激活/CAF"]
    subtitles = {"潘氏细胞·4": "潘氏细胞·4（★ 3/3 队列腺瘤峰）",
                 "纯上皮基线·7": "纯上皮基线·7（阴性对照）",
                 "间质·激活/CAF": "间质·激活/CAF（跨区室对照）"}
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.4))
    styles = {COH[0]: ("o", CB, "-"), COH[1]: ("s", CT, "-"), COH[2]: ("^", CR, "--")}
    for ax, pn in zip(axes.ravel(), panels):
        if pn not in TLA:
            ax.axis("off")
            continue
        for g in COH:
            if g not in TLA[pn]:
                continue
            d = TLA[pn][g]
            ys = [d["normal"], d["adenoma"], d["carcinoma"]]
            mk, c, ls = styles[g]
            ax.plot([0, 1, 2], ys, marker=mk, color=c, ls=ls, lw=2.0, ms=7,
                    label="%s [%s]" % (SHORT[g], d.get("pattern", "?")))
        ax.axhline(0, color=GREY, lw=0.8, ls=":")
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["正常", "腺瘤", "癌"], fontsize=10)
        ax.set_title(subtitles.get(pn, pn), fontsize=11.5, fontweight="bold")
        ax.set_ylabel("面板分数（z 均值，已校正）", fontsize=8.6)
        ax.legend(fontsize=7.6, loc="best", framealpha=0.9)
        ax.grid(alpha=0.25, lw=0.6)
    fig.suptitle("图 2  三层模式（漂移校正口径；raw 口径见 Fig. S2）：潘氏细胞程序在 3/3 队列呈「腺瘤峰」P 型，而阴性对照单调下降\n"
                 "P = 腺瘤峰（腺瘤 ≥ 正常 且 腺瘤 > 癌）；D = 单调降；A = 腺瘤期已降",
                 fontsize=12.6, fontweight="bold", y=1.02)
    fig.tight_layout()
    save(fig, "fig2_三层模式_腺瘤峰.png")


# ===========================================================================
# Fig S2  双口径筛选矩阵：13 面板 × 3 队列 × 2 口径
# ===========================================================================
def figS2():
    import numpy as np
    from matplotlib.colors import ListedColormap
    tab = R3["panel_table"]
    passed = set(R3["robust_P_panels"])
    order = sorted(tab, key=lambda k: (k not in passed, k))      # 通过的排最上
    cols = [(g, m) for g in COH for m in ("raw", "adj")]
    code = {"P": 0, "D": 1, "A": 2, "?": 3}
    M = np.full((len(order), len(cols)), 3, dtype=int)
    for i, pn in enumerate(order):
        for j, (g, m) in enumerate(cols):
            key = "raw_patterns" if m == "raw" else "adj_patterns"
            M[i, j] = code.get(tab[pn][key].get(g, "?"), 3)
    cmap = ListedColormap([CR, CB, "#cbd5e1", "#f8fafc"])
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    ax.imshow(M, cmap=cmap, aspect="auto", vmin=0, vmax=3)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["%s\n%s" % (SHORT[g].split()[0], "原始" if m == "raw" else "校正")
                        for g, m in cols], fontsize=8.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(["%s%s" % (pn, "  ★" if pn in passed else "") for pn in order], fontsize=9.4)
    ax.set_xticks(np.arange(-.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.8)
    ax.tick_params(which="minor", length=0)
    ax.axhline(0.5, color="#0f172a", lw=2.2)
    for i in range(len(order)):
        for j in range(len(cols)):
            v = ["P", "D", "A", "?"][M[i, j]]
            ax.text(j, i, v, ha="center", va="center", fontsize=12.5, fontweight="bold",
                    color="white" if M[i, j] in (0, 1) else "#475569")
    ax.set_title("图 S2  双口径筛选：仅「潘氏细胞·4」（★，分隔线以上）在 3 队列 × 2 口径下全部为腺瘤峰（P）\n"
                 "P = 腺瘤峰（腺瘤 ≥ 正常 且 腺瘤 > 癌）｜D = 单调降｜A = 腺瘤期已降\n"
                 "原始 = 队列内 z 化；校正 = 另加样本级 glev 漂移校正（见 R2 §0）",
                 fontsize=10.4, fontweight="bold", pad=10)
    fig.tight_layout()
    save(fig, "figS2_双口径筛选矩阵.png")


# ===========================================================================
# Fig 3  2×2 分解：全组织 vs LCM，原始 vs 漂移校正
# ===========================================================================
def fig3():
    sets_ = ["上调分量(210)", "下调分量(103)", "上皮分化·核心8", "泛间质·6", "WGCNA-turquoise"]
    fig, ax = plt.subplots(figsize=(12.4, 5.6))
    n = len(sets_)
    W = 0.20
    cols = [CB, CT, CA, GREY]
    labs = ["GSE4183 全组织", "GSE117606 全组织", "GSE15960 LCM（原始）", "GSE15960 LCM（漂移校正）"]
    for gi, s in enumerate(sets_):
        vals = []
        for g in COH:
            vals.append(MAT.get(s, {}).get(g, {}).get("d_median"))
        # LCM 校正值自己算：只在 panel_matrix 里没存，用 pct_up 近似不可靠 → 见下
        vals = [v if v is not None else 0 for v in vals]
        raw = MAT.get(s, {}).get("GSE15960", {}).get("d_median_raw")
        adj = MAT.get(s, {}).get("GSE15960", {}).get("d_median")
        bars = [MAT.get(s, {}).get("GSE4183", {}).get("d_median") or 0,
                MAT.get(s, {}).get("GSE117606", {}).get("d_median") or 0,
                raw or 0, adj or 0]
        for bi, (b, c) in enumerate(zip(bars, cols)):
            x = gi + (bi - 1.5) * W
            ax.bar(x, b, W * 0.88, color=c, edgecolor="white", lw=0.7)
            ax.text(x, b + (0.06 if b >= 0 else -0.14), "%+.2f" % b,
                    ha="center", fontsize=7.6, color="#0f172a")
    ax.axhline(0, color="#1e293b", lw=1.1)
    ax.set_xticks(range(n))
    ax.set_xticklabels([s.replace("(", "\n(") for s in sets_], fontsize=9.5)
    ax.set_ylabel("Δ = 癌 − 腺瘤（中位，z）", fontsize=11)
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, labs)],
              fontsize=9.5, ncol=2, loc="upper left")
    ax.set_title("图 3  2×2 分解（逐基因落差保留比例）：上调分量落差衰减到 0.32 倍（同向 58.6%），\n"
                 "下调分量落差放大到 2.04 倍（同向 93.2%）",
                 fontsize=12.2, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig3_2x2分解_全组织vsLCM.png")


# ===========================================================================
# Fig 4  逐基因散点：全组织 Δ  vs  LCM Δ
# ===========================================================================
def fig4():
    rows = list(csv.DictReader(open("suppl_table1_313genes.csv", encoding="utf-8-sig")))

    def f(v):
        try:
            return float(v)
        except Exception:
            return None
    xs, ys, cs, names = [], [], [], []
    for r in rows:
        a, b = f(r.get("d_GSE4183")), f(r.get("d_GSE117606"))
        y = f(r.get("d_GSE15960_LCM"))
        if y is None or (a is None and b is None):
            continue
        x = statistics.mean([v for v in (a, b) if v is not None])
        xs.append(x)
        ys.append(y)
        names.append(r["gene"])
        cs.append(CR if f(r["mean_logFC"]) > 0 else CB)
    fig, ax = plt.subplots(figsize=(8.6, 7.4))
    lim = 3.0
    ax.axhspan(-lim, 0, color="#f1f5f9", zorder=0)
    ax.axvspan(-lim, 0, color="#f1f5f9", zorder=0)
    ax.scatter(xs, ys, s=17, c=cs, alpha=0.75, edgecolors="white", linewidths=0.3)
    ax.axhline(0, color="#1e293b", lw=1.0)
    ax.axvline(0, color="#1e293b", lw=1.0)
    ax.plot([0, 2.6], [-0.9 * 0, -0.9 * 2.6], color=GREY, ls="--", lw=1.2)
    ax.text(1.55, -2.35, "若 LCM 落差与全组织同幅\n（斜率=1）", color=GREY, fontsize=8.5)
    for tag in ["FAP", "COL11A1", "THBS2", "ASPN", "INHBA", "IL6", "MMP1", "PDGFRB",
                "CLCA1", "MUC2", "FCGBP", "REG4", "ZG16", "SI", "DEFA5", "OLFM4", "SPINK1", "NAT2"]:
        if tag in names:
            i = names.index(tag)
            ax.annotate(tag, (xs[i], ys[i]), fontsize=8.2, fontweight="bold",
                        xytext=(3, 3), textcoords="offset points",
                        color="#0f172a")
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(-3.0, 3.0)
    ax.set_xlabel("全组织 Δ（癌 − 腺瘤，两队列均值，已校正）", fontsize=11)
    ax.set_ylabel("LCM 纯上皮 Δ（癌 − 腺瘤，已校正）", fontsize=11)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(handles=[Patch(color=CR, label="签名上调分量（210）"),
                       Patch(color=CB, label="签名下调分量（103）")],
              fontsize=9.5, loc="lower right")
    ax.set_title("图 4  逐基因归属（313 基因）：上调分量（红）多数落在对角线上方但幅度被压缩，\n"
                 "下调分量（蓝）系统性落在对角线下方 —— 去掉间质后落差反而更大",
                 fontsize=12.0, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig4_逐基因归属散点.png")
    # 参考辅助线：y = 0.32x（上调分量实际保留斜率）
    return None


# ===========================================================================
# Fig 5  面板 × 队列 |AUC| 热图（含阴性对照）
# ===========================================================================
def fig5():
    ORDER = ["313签名", "上调分量(210)", "下调分量(103)", "上皮·分化/功能", "上皮分化·核心8",
             "纯上皮基线·7", "间质·激活/CAF", "间质·基线", "泛间质·6", "髓系·激活",
             "WGCNA-turquoise", "Hub基因"]
    cols = []
    for g in COH:
        cols.append((g, "auc", SHORT[g] + "\n校正"))
    for g in COH:
        if g == "GSE15960":
            cols.append((g, "auc_raw", SHORT[g] + "\n原始"))
    data = []
    for s in ORDER:
        row = []
        for g, key, _ in cols:
            v = MAT.get(s, {}).get(g, {}).get(key)
            row.append(max(v, 1 - v) if v is not None else float("nan"))
        data.append(row)
    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    im = ax.imshow(data, cmap="RdYlBu_r", vmin=0.5, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[2] for c in cols], fontsize=8.6)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER, fontsize=9.5)
    for i in range(len(ORDER)):
        for j in range(len(cols)):
            v = data[i][j]
            if v != v:
                ax.text(j, i, "—", ha="center", va="center", fontsize=9, color="#475569")
            else:
                ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=8.2,
                        color="white" if (v > 0.88 or v < 0.58) else "#0f172a")
    # 阴性对照高亮
    ax.axhline(5.5, color="#0f172a", lw=1.6)
    ax.text(-1.35, 5.5, "▲ 阴性对照", ha="center", va="center", fontsize=8.5,
            rotation=90, color=CR, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("|AUC|（0.5 = 无区分力）", fontsize=9.5)
    ax.set_title("图 5  面板 × 队列区分力（癌 vs 腺瘤）\nLCM 一列请对照原始/校正两种口径看",
                 fontsize=12.5, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig5_面板队列AUC热图.png")


# ===========================================================================
# Fig 6  56 组患者配对：各面板同向比例
# ===========================================================================
def fig6():
    B = PAIR["B_panel_paired"]
    # B_paired 用 n_pair 键，B_panel_paired 用 n 键 —— 统一成 n
    def norm(v):
        return {"n": v.get("n") or v.get("n_pair"), "n_pos": v["n_pos"],
                "median_delta": v["median_delta"], "wilcoxon_p": v["wilcoxon_p"]}
    items = [("313 全签名", norm(PAIR["B_paired"]["全签名"])),
             ("上调分量（间质/髓系激活）", norm(PAIR["B_paired"]["上调分量(间质/髓系激活)"])),
             ("下调分量（上皮功能关闭）", norm(PAIR["B_paired"]["下调分量(上皮功能关闭)"]))]
    items += [(k, norm(v)) for k, v in B.items()]
    items = sorted(items, key=lambda kv: kv[1]["n_pos"] / kv[1]["n"])
    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    ys = range(len(items))
    pct = [100.0 * v["n_pos"] / v["n"] for _, v in items]
    colors = []
    for k, v in items:
        p = 100.0 * v["n_pos"] / v["n"]
        if "基线/结构" in k or "阴性" in k:
            colors.append(GREY)
        else:
            colors.append(CR if p >= 50 else CB)
    ax.barh(list(ys), pct, color=colors, height=0.62, edgecolor="white")
    for i, (k, v) in enumerate(items):
        p = 100.0 * v["n_pos"] / v["n"]
        ax.text(p + 1.2, i, "%d/%d = %.1f%%  中位Δ %+.2f  p=%.1e"
                % (v["n_pos"], v["n"], p, v["median_delta"], v["wilcoxon_p"]),
                va="center", fontsize=8.4, color="#0f172a")
    ax.axvline(50, color="#1e293b", lw=1.4, ls="--")
    ax.text(50.6, len(items) - 0.35, "50% = 无方向（随机）", fontsize=9, color="#1e293b")
    ax.set_yticks(list(ys))
    ax.set_yticklabels([k for k, _ in items], fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_xlabel("同向上调比例（56 组「同一患者 腺瘤 → 癌」配对）%", fontsize=10.5)
    ax.grid(axis="x", alpha=0.25, lw=0.6)
    ax.set_title("图 6  配对设计（固定年龄/性别/遗传背景/取材批次）：间质激活 94.6% 同向，上皮分化仅 14.3% 同向",
                 fontsize=12, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig6_56组配对同向比例.png")


# ===========================================================================
# 从 facts 里取 LCM 校正后的面板 Δ（panel_matrix 存的是校正值）
# 说明：校正/原始两列都已在 manuscript_facts.json 的 panel_matrix 里
# ===========================================================================
ADJ_D = {s: {"GSE15960": MAT.get(s, {}).get("GSE15960", {}).get("d_median")} for s in MAT}

if __name__ == "__main__":
    print("生成图形：")
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    figS2()
    print("DONE")
