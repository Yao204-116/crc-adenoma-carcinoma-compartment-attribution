# -*- coding: utf-8 -*-
"""
发表级图形生成（matplotlib）
================================================================================
数据源（全部为本轮实测）：
  three_layer_pattern.json      三层（正常/腺瘤/癌）面板分数 + 模式标签
  manuscript_facts.json         面板×队列 AUC（原始 & 漂移校正）、全局漂移
  suppl_table1_313genes.csv     313 基因逐基因 Δ（三队列，已校正）
  gse117606_paired_results.json 56 组患者配对

输出 figs/ 下 6 张 300dpi PNG：Figure 2 / 3 / 4 / 5 / 6 / S2
用法：python make_figs.py

★ 英文化铁律（见 fig_en.py 顶部说明）：**数据键一律保持中文**
  （`panels` / `ORDER` / `sets_` / `order` 既查字典又上屏），
  只在 `set_title` / `set_xticklabels` / `set_yticklabels` 这类**文字上屏行**上过
  `fig_en.panel_en()`。保存统一走 `fig_en.finish()`，它会在落盘前扫一遍画布，
  残留汉字就报错退出。
"""
import csv, json, os, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import fig_en

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)

fig_en.use_english_font()
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
# 显示用短名（值只用于上屏，键仍是队列代码）—— 用词对齐稿件图注
SHORT = {"GSE4183": "GSE4183 whole tissue", "GSE117606": "GSE117606 whole tissue",
         "GSE15960": "GSE15960 LCM"}


def save(fig, name):
    fig_en.finish(fig, os.path.join(FIG, name), bbox_inches="tight", facecolor="white")


# ===========================================================================
# Figure 6  三层模式：腺瘤峰（6 面板 × 3 队列）+ 阴性对照
# ===========================================================================
def fig2():
    # ★ 用漂移校正口径（保守口径）；FACTS["three_layer"]["raw"] 为原始口径
    TLA = FACTS["three_layer"]["adj"]
    # 键名必须与数据完全一致（中文）；上屏时才翻英文
    panels = ["潘氏细胞·4", "黏液杯状·5", "上皮分化·核心8", "相代谢·5",
              "纯上皮基线·7", "间质·激活/CAF"]
    subtitles = {"潘氏细胞·4": "Paneth-cell·4\n(★ adenoma peak in 3/3 cohorts)",
                 "纯上皮基线·7": "pure-epithelial-baseline·7\n(negative control)",
                 "间质·激活/CAF": "stroma/activated-CAF\n(cross-compartment control)"}
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
        ax.set_xticklabels(["Normal", "Adenoma", "Carcinoma"], fontsize=10)
        ax.set_title(fig_en.panel_en(subtitles.get(pn, pn)), fontsize=11.5, fontweight="bold")
        ax.set_ylabel("Panel score (mean z, drift-corrected)", fontsize=8.6)
        ax.legend(fontsize=7.6, loc="best", framealpha=0.9)
        ax.grid(alpha=0.25, lw=0.6)
    fig.suptitle(fig_en.wrap(
        "Figure 6  Three-layer pattern (drift-corrected caliber; raw caliber in Fig. S2): "
        "the Paneth-cell programme shows an adenoma peak (P) in 3/3 cohorts, whereas the "
        "negative control declines monotonically\n"
        "P = adenoma peak (adenoma ≥ normal and adenoma > carcinoma) · D = monotonic decrease · "
        "A = already decreased at adenoma", 12.8, 12.0),
        fontsize=12.0, fontweight="bold", y=1.02)
    fig.tight_layout()
    save(fig, "fig6_three_layer_pattern_adenoma_peak.png")


# ===========================================================================
# Figure S2  双口径筛选矩阵：13 面板 × 3 队列 × 2 口径
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
    ax.set_xticklabels(["%s\n%s" % (SHORT[g].split()[0], "Raw" if m == "raw" else "Corrected")
                        for g, m in cols], fontsize=8.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(["%s%s" % (fig_en.panel_en(pn), "  ★" if pn in passed else "")
                        for pn in order], fontsize=9.4)
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
    ax.set_title(fig_en.wrap(
        "Figure S2  Dual-caliber screening matrix (13 panels × 3 cohorts × 2 calibers): "
        "only Paneth-cell·4 (★, above the divider) is an adenoma peak (P) in all three cohorts "
        "and in both calibers\n"
        "P = adenoma peak (adenoma ≥ normal and adenoma > carcinoma) | D = monotonic decrease | "
        "A = already decreased at adenoma\n"
        "Raw = within-cohort z; Corrected = plus sample-level glev drift correction", 6.0, 10.0),
        fontsize=10.0, fontweight="bold", pad=10)
    fig.tight_layout()
    save(fig, "figS2_double_caliber_screening_matrix.png")


# ===========================================================================
# Figure 3  2×2 分解：全组织 vs LCM，原始 vs 漂移校正
# ===========================================================================
def fig3():
    sets_ = ["上调分量(210)", "下调分量(103)", "上皮分化·核心8", "泛间质·6", "WGCNA-turquoise"]
    fig, ax = plt.subplots(figsize=(12.4, 5.6))
    n = len(sets_)
    W = 0.20
    cols = [CB, CT, CA, GREY]
    labs = ["GSE4183 whole tissue", "GSE117606 whole tissue",
            "GSE15960 LCM (raw)", "GSE15960 LCM (drift-corrected)"]
    for gi, s in enumerate(sets_):
        # LCM 那一列直接读 panel_matrix：d_median_raw 是原始口径，d_median 是校正口径
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
    ax.set_xticklabels([fig_en.panel_en(s) for s in sets_], fontsize=9.5)
    ax.set_ylabel("Δ = carcinoma − adenoma (median, z)", fontsize=11)
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, labs)],
              fontsize=9.5, ncol=2, loc="upper left")
    ax.set_title(fig_en.wrap(
        "Figure 3  2×2 decomposition (per-gene retention of the carcinoma−adenoma gap): "
        "the up-component decays to 0.32× (58.6% concordant), the down-component is amplified "
        "2.04× (93.2% concordant)", 11.2, 11.5),
        fontsize=11.5, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig3_2x2_decomposition_whole_tissue_vs_LCM.png")


# ===========================================================================
# Figure 4  逐基因散点：全组织 Δ  vs  LCM Δ
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
    ax.text(-2.85, -2.72, "grey dashed line: y = −0.9x\n(if the LCM gap equalled whole tissue)",
            ha="left", va="bottom", color=GREY, fontsize=8.5)
    anns = []       # 基因名标注：稍后交给 fig_en.declutter 自动错开
    for tag in ["FAP", "COL11A1", "THBS2", "ASPN", "INHBA", "IL6", "MMP1", "PDGFRB",
                "CLCA1", "MUC2", "FCGBP", "REG4", "ZG16", "SI", "DEFA5", "OLFM4", "SPINK1", "NAT2"]:
        if tag in names:
            i = names.index(tag)
            anns.append(ax.annotate(tag, (xs[i], ys[i]), fontsize=8.2, fontweight="bold",
                                    xytext=(3, 3), textcoords="offset points",
                                    color="#0f172a"))
    ax.set_xlim(-3.0, 3.0)
    ax.set_ylim(-3.0, 3.0)
    ax.set_xlabel("whole-tissue Δ (carcinoma − adenoma, mean of two cohorts, drift-corrected)", fontsize=11)
    ax.set_ylabel("LCM pure-epithelium Δ (carcinoma − adenoma, drift-corrected)", fontsize=11)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(handles=[Patch(color=CR, label="signature up-component (210)"),
                       Patch(color=CB, label="signature down-component (103)")],
              fontsize=9.5, loc="lower right")
    ax.set_title(fig_en.wrap(
        "Figure 4  Per-gene attribution (313 genes): most up-component genes (red) lie above the "
        "diagonal but with compressed amplitude, down-component genes (blue) lie systematically "
        "below it — removing the stroma makes the gap larger", 7.2, 10.5),
        fontsize=10.5, fontweight="bold", pad=12)
    fig.tight_layout()
    # ★ 英文基因名仍会互压（FAP×COL11A1 81%、CLCA1×FCGBP 39%）→ 自动沿 y 错开
    fig_en.declutter(fig, anns)
    save(fig, "fig4_per_gene_attribution_scatter.png")
    # 参考辅助线：y = 0.32x（上调分量实际保留斜率）
    return None


# ===========================================================================
# Figure 5  面板 × 队列 |AUC| 热图（含阴性对照）
# ===========================================================================
def fig5():
    ORDER = ["313签名", "上调分量(210)", "下调分量(103)", "上皮·分化/功能", "上皮分化·核心8",
             "纯上皮基线·7", "间质·激活/CAF", "间质·基线", "泛间质·6", "髓系·激活",
             "WGCNA-turquoise", "Hub基因"]
    cols = []
    for g in COH:
        cols.append((g, "auc", SHORT[g] + "\nCorrected"))
    for g in COH:
        if g == "GSE15960":
            cols.append((g, "auc_raw", SHORT[g] + "\nRaw"))
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
    ax.set_yticklabels([fig_en.panel_en(s) for s in ORDER], fontsize=9.5)
    # ★ 阴性对照行：染红加粗 + 标题里点明。
    #   原来是在轴外画一条竖排的「▲ negative control」——中文面板名换英文后变长，
    #   左侧已无空白，那段竖排文字会直接压在 y 轴标签上（目视复查发现）。
    for t, s in zip(ax.get_yticklabels(), ORDER):
        if s == "纯上皮基线·7":
            t.set_color(CR)
            t.set_fontweight("bold")
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
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("|AUC| (0.5 = no discrimination)", fontsize=9.5)
    ax.set_title("Figure 5  Panel × cohort discrimination (carcinoma vs adenoma)\n"
                 "For the LCM column, read the raw and drift-corrected calibers side by side\n"
                 "red label = negative control (pure-epithelial-baseline·7)",
                 fontsize=12.5, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig5_panel_cohort_AUC_heatmap.png")


# ===========================================================================
# Figure 2  56 组患者配对：各面板同向比例
# ===========================================================================
def fig6():
    B = PAIR["B_panel_paired"]
    # B_paired 用 n_pair 键，B_panel_paired 用 n 键 —— 统一成 n
    def norm(v):
        return {"n": v.get("n") or v.get("n_pair"), "n_pos": v["n_pos"],
                "median_delta": v["median_delta"], "wilcoxon_p": v["wilcoxon_p"]}
    items = [("313-gene signature", norm(PAIR["B_paired"]["全签名"])),
             ("stromal/myeloid-activation component", norm(PAIR["B_paired"]["上调分量(间质/髓系激活)"])),
             ("epithelial-function component", norm(PAIR["B_paired"]["下调分量(上皮功能关闭)"]))]
    items += [(k, norm(v)) for k, v in B.items()]
    items = sorted(items, key=lambda kv: kv[1]["n_pos"] / kv[1]["n"])
    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    ys = range(len(items))
    pct = [100.0 * v["n_pos"] / v["n"] for _, v in items]
    colors = []
    for k, v in items:
        p = 100.0 * v["n_pos"] / v["n"]
        # ★ k 仍是中文数据键：这里的子串判断不能改，否则阴性对照的灰底会丢
        if "基线/结构" in k or "阴性" in k:
            colors.append(GREY)
        else:
            colors.append(CR if p >= 50 else CB)
    ax.barh(list(ys), pct, color=colors, height=0.62, edgecolor="white")
    for i, (k, v) in enumerate(items):
        p = 100.0 * v["n_pos"] / v["n"]
        ax.text(p + 1.2, i, "%d/%d = %.1f%%  median Δ %+.2f  p=%.1e"
                % (v["n_pos"], v["n"], p, v["median_delta"], v["wilcoxon_p"]),
                va="center", fontsize=8.4, color="#0f172a")
    ax.axvline(50, color="#1e293b", lw=1.4, ls="--")
    ax.text(50.6, len(items) - 0.35, "50% = no direction (random)", fontsize=9, color="#1e293b")
    ax.set_yticks(list(ys))
    ax.set_yticklabels([fig_en.panel_en(k) for k, _ in items], fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Same-direction proportion (%) across 56 matched patient pairs (adenoma → carcinoma)", fontsize=10.5)
    ax.grid(axis="x", alpha=0.25, lw=0.6)
    ax.set_title("Figure 2  Paired design (age, sex, germline background and sampling batch held fixed):\n"
                 "stromal activation 94.6% concordant, epithelial differentiation only 14.3% concordant",
                 fontsize=12, fontweight="bold", pad=12)
    fig.tight_layout()
    save(fig, "fig2_56_pair_concordance.png")


if __name__ == "__main__":
    print("Writing figures:")
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()
    figS2()
    print("DONE")
