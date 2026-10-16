# -*- coding: utf-8 -*-
"""
出路 B 第一步：跨队列共识签名（真正跑，不是方案）
4 个独立队列，同一生物学对比方向（腺瘤/息肉 → 癌）：
  GSE41258   Tumor vs Polyp       GPL96/97     n=234
  GSE37364   Carcinoma vs Adenoma GPL570
  GSE20916   Carcinoma vs Adenoma GPL570 (macro)
  GSE41657   Carcinoma vs Adenoma Agilent
方法：
  ① 基因级方向一致性（sign concordance）
  ② Stouffer 加权合并 z（用各队列 p 值与 logFC 方向）
  ③ Fisher 法合并 p
  ④ 分层：Core（4/4 同向显著）/ High（≥3/4 同向 + meta FDR<0.05）
输出：共识基因表 + 与 WGCNA/富集结果的交集 + 细胞来源归属
"""
import csv, json, math, os, statistics
from collections import Counter

# 定位策略：默认取脚本所在目录，可用 AACA_BASE 覆盖，便于从任意 cwd 调用。
BASE = os.environ.get("AACA_BASE") or os.path.dirname(os.path.abspath(__file__))
# 上游 R 流程产物（limma DEG / WGCNA / 富集交集）。归档副本随仓库发布在 upstream_r/ 下；
# 该目录不存在时（本机开发环境）回退到原始位置 D:/2。
_UPSTREAM = os.environ.get("AACA_UPSTREAM_R") or os.path.join(BASE, "upstream_r")
D2 = _UPSTREAM if os.path.isdir(_UPSTREAM) else r"D:/2"

COHORTS = [
    ("GSE41258", "Tumor vs Polyp", os.path.join(D2, "01.limma/DEG_all_Tumor_vs_Polyp.csv"), "name"),
    ("GSE37364", "Carcinoma vs Adenoma",
     os.path.join(BASE, "02.result/GSE37364_DEG_Carcinoma_vs_Adenoma.csv"), "symbol"),
    ("GSE20916", "Carcinoma vs Adenoma (macro)",
     os.path.join(BASE, "02.result/GSE20916macro_DEG_Carcinoma_vs_Adenoma.csv"), "symbol"),
    ("GSE41657", "Carcinoma vs Adenoma (Agilent)",
     os.path.join(BASE, "02.result/GSE41657_DEG_Carcinoma_vs_Adenoma.csv"), "symbol"),
]

# 细胞类型面板（沿用上一轮，用于把共识基因落到细胞来源）
PANELS = {
    "上皮·基线/结构": "EPCAM KRT8 KRT18 KRT19 CDH1 CLDN3 CLDN4 CLDN7 DSP JUP PKP3 ITGB4 ITGA6 SFN RAB25 PRSS8 GRHL2 ELF3 EHF ESRP1 ESRP2 DSC2 PKP2 OCLN ANXA1 ANXA2 S100A11 GSN",
    "上皮·分化/功能": "MUC2 MUC4 MUC12 MUC13 MUC20 TFF1 TFF2 TFF3 SI OLFM4 PIGR ZG16 CA1 CA2 CEACAM7 GUCA2A GUCA2B AQP8 SLC26A3 SLC51B HPGD CD177 CLCA4 PLAC8 FABP1 KRT20 MEP1A ANPEP APOA4 CDX2 VIL1 REG4 SPINK4 FCGBP AGR2 ITLN1 HNF4A KRT23 ADH1C SLC12A2 MUC1",
    "间质·基线": "DCN LUM PDGFRA FBLN1 FBLN2 OGN PRELP FMOD C7 CFH SERPINF1 PCOLCE MFAP4 MFAP5 ECM1 PTGDS SOD3 MMP2 IGFBP7 CCN1 CCN2 IGFBP5 GAS6 LAMA4 LAMB1 ITGA1 CAV1 CAV2 RARRES2 NNMT NBL1",
    "间质·激活/CAF": "FAP THBS2 COL10A1 COL11A1 INHBA COMP TNC LRRC15 CTHRC1 SULF1 SULF2 GREM1 WNT5A POSTN ASPN MMP11 MMP13 MMP1 LOXL2 CDH11 PXDN ADAM12 SEMA3C NID2 ADAMTS2 ISLR EMILIN1 COL8A1 SPARC THY1 SERPINE1 TIMP1 MMP3 MMP14 LOX PLOD2 BGN FNDC1 COL16A1 THBS1 CTGF ANGPTL4 VEGFA",
    "免疫·T/B 基线": "CD3D CD3E CD3G CD2 CD247 TRAC TRBC2 CD8A CD8B CD4 IL7R CD28 LCK ZAP70 THEMIS CCR7 SELL LTB MS4A1 CD79A CD79B CD19 BANK1 IGHD TNFRSF13B SKAP1 ITK TRAT1 CD27 CD52 PTPRC",
    "免疫·激活/效应": "GZMA GZMB GZMK GZMH PRF1 NKG7 GNLY KLRD1 CTLA4 PDCD1 LAG3 TIGIT FOXP3 IL2RA CD69 ICOS TNFRSF4 TNFRSF18 CD40LG IFNG CCL5 CXCL9 CXCL10 CXCL13 IGHG1 IGHG2 IGHG3 IGHG4 IGKC JCHAIN MZB1 SDC1",
    "髓系·基线": "CD68 CD163 MRC1 CSF1R ITGAM FCGR2A CD14 LYZ TYROBP AIF1 FCER1G CTSB CTSS CTSD MAF MAFB MS4A7 VSIG4 FCGR1A SLC11A1",
    "髓系·激活": "SPP1 CCL18 CCL2 CCL3 CCL4 CCL20 APOE TREM2 MMP9 MMP12 SIGLEC1 CXCL8 IL10 IDO1 CD274 PDCD1LG2 LGALS3 S100A8 S100A9 S100A12 FCN1 C5AR1 MNDA FCGR3B MMP8 TCN1 CCL19",
}
PANELS = {k: set(v.split()) for k, v in PANELS.items()}


def load(path, symcol):
    """→ {gene: (logFC, p)}  多探针取 |logFC| 最大者"""
    if not os.path.exists(path):
        return None
    rows = list(csv.DictReader(open(path, encoding="utf-8-sig")))
    if not rows:
        return None
    cols = list(rows[0].keys())
    if symcol not in cols:
        for c in ("symbol", "name", "gene", "", "X"):
            if c in cols:
                symcol = c
                break
        else:
            symcol = cols[0]
    out = {}
    for r in rows:
        s = (r.get(symcol) or "").strip().strip('"')
        if not s or s in ("NA", "."):
            continue
        try:
            fc = float(r["logFC"])
        except (KeyError, ValueError, TypeError):
            continue
        p = 1.0
        for k in ("P.Value", "p", "Pvalue"):
            if k in r and r[k] not in (None, ""):
                try:
                    p = float(r[k])
                    break
                except ValueError:
                    pass
        p = min(max(p, 1e-300), 1.0)
        if s not in out or abs(fc) > abs(out[s][0]):
            out[s] = (fc, p)
    return out


def phi_inv(q):
    return statistics.NormalDist().inv_cdf(min(max(q, 1e-300), 1 - 1e-16))


def fdr(ps):
    """Benjamini-Hochberg"""
    n = len(ps)
    idx = sorted(range(n), key=lambda i: ps[i])
    out = [1.0] * n
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = idx[rank]
        val = ps[i] * n / (rank + 1)
        prev = min(prev, val)
        out[i] = min(prev, 1.0)
    return out


print("=" * 112)
print("出路 B 第一步：跨 4 队列腺瘤/息肉 → 癌 共识签名")
print("=" * 112)

data = {}
for name, contrast, path, symcol in COHORTS:
    d = load(path, symcol)
    if d is None:
        print("!! 缺失:", name, path)
        continue
    data[name] = d
    nsig = sum(1 for fc, p in d.values() if p < 0.05)
    print("  %-10s %-30s 基因 %6d   其中 p<0.05 %5d" % (name, contrast, len(d), nsig))

names = [c[0] for c in COHORTS if c[0] in data]
print("\n参与合并的队列: %s (%d 个)" % (", ".join(names), len(names)))

# ---- 队列间 logFC 相关性（QC）----
print("\n" + "-" * 112)
print("队列间一致性 QC（共同基因的 logFC Pearson / Spearman）")
print("-" * 112)


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / (sx * sy)


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = data[names[i]], data[names[j]]
        common = sorted(set(a) & set(b))
        xs = [a[g][0] for g in common]
        ys = [b[g][0] for g in common]
        print("  %-10s ~ %-10s  n=%6d   Pearson=%+.3f   Spearman=%+.3f"
              % (names[i], names[j], len(common), pearson(xs, ys), spearman(xs, ys)))

# ---- 合并 ----
rows = []
allgenes = set()
for n in names:
    allgenes |= set(data[n])

for g in allgenes:
    present = [(n, data[n][g]) for n in names if g in data[n]]
    k = len(present)
    if k < 2:
        continue
    fcs = [v[0] for _, v in present]
    ps = [v[1] for _, v in present]
    mean_fc = statistics.mean(fcs)
    med_fc = statistics.median(fcs)
    # 方向一致性
    dirs = [1 if f > 0 else (-1 if f < 0 else 0) for f in fcs]
    mdir = 1 if med_fc > 0 else -1
    agree = sum(1 for d in dirs if d == mdir)
    # Stouffer
    zi = [phi_inv(1 - p / 2.0) * (1 if f > 0 else -1) for f, p in zip(fcs, ps)]
    z_comb = sum(zi) / math.sqrt(k)
    p_stouffer = math.erfc(abs(z_comb) / math.sqrt(2))
    # Fisher
    chi = -2.0 * sum(math.log(p) for p in ps)
    # 卡方(2k) 上尾 —— 用 Wilson-Hilferty 近似，避免依赖 scipy
    df = 2 * k
    zc = ((chi / df) ** (1.0 / 3.0) - (1 - 2.0 / (9 * df))) / math.sqrt(2.0 / (9 * df))
    p_fisher = 0.5 * math.erfc(zc / math.sqrt(2))
    min_p = min(ps)
    rows.append({"gene": g, "n_cohort": k, "mean_logFC": mean_fc, "median_logFC": med_fc,
                 "agree": agree, "dir_ratio": agree / float(k), "z_stouffer": z_comb,
                 "p_stouffer": p_stouffer, "chi2": chi, "p_fisher": p_fisher, "min_p": min_p,
                 "fcs": [round(f, 3) for f in fcs]})

print("\n可合并基因总数（出现在 ≥2 个队列）: %d" % len(rows))

# meta FDR（对 Stouffer p 做 BH）
ps_all = [r["p_stouffer"] for r in rows]
qs = fdr(ps_all)
for r, q in zip(rows, qs):
    r["q_stouffer"] = q

# ---- 分层 ----
core = [r for r in rows if r["n_cohort"] == len(names) and r["agree"] == r["n_cohort"]
        and r["min_p"] < 0.05 and r["q_stouffer"] < 0.05]
high = [r for r in rows if r["n_cohort"] >= max(3, len(names) - 1) and r["dir_ratio"] == 1.0
        and r["q_stouffer"] < 0.05]
core.sort(key=lambda r: -abs(r["mean_logFC"]))
high.sort(key=lambda r: -abs(r["mean_logFC"]))

print("\n" + "=" * 112)
print("共识分层结果")
print("=" * 112)
print("  Core 共识（%d/4 队列全同向 + 每队列 p<0.05 + meta FDR<0.05）: %d 个基因" % (len(names), len(core)))
print("  High 共识（≥3 队列全同向 + meta FDR<0.05）            : %d 个基因" % len(high))
up_core = sum(1 for r in core if r["mean_logFC"] > 0)
print("     其中 Core 上调 %d / 下调 %d" % (up_core, len(core) - up_core))

print("\nCore 共识 Top 30（按 |mean logFC| 排序）：")
print("  %-12s %6s %8s %8s %6s  %s" % ("gene", "k", "meanFC", "medFC", "agree", "各队列 logFC"))
for r in core[:30]:
    print("  %-12s %6d %+8.3f %+8.3f %6s  %s"
          % (r["gene"], r["n_cohort"], r["mean_logFC"], r["median_logFC"],
             "%d/%d" % (r["agree"], r["n_cohort"]), r["fcs"]))

# ---- 细胞来源归属（排他性）----
print("\n" + "=" * 112)
print("Core 共识基因的细胞来源归属（排他性：只被单一面板命中）")
print("=" * 112)
excl = {k: [] for k in PANELS}
shared = []
unassigned = []
for r in core:
    g = r["gene"]
    owners = [k for k, s in PANELS.items() if g in s]
    if len(owners) == 1:
        excl[owners[0]].append(r)
    elif len(owners) > 1:
        shared.append((g, owners))
    else:
        unassigned.append(r)

for k in PANELS:
    lst = excl[k]
    if not lst:
        print("  %-16s 0" % k)
        continue
    ups = [r for r in lst if r["mean_logFC"] > 0]
    mfc = statistics.median([r["mean_logFC"] for r in lst])
    print("  %-16s %3d 命中 | 上调 %2d / 下调 %2d | 中位 meanFC %+6.3f | %s"
          % (k, len(lst), len(ups), len(lst) - len(ups), mfc,
             ", ".join("%s(%+.2f)" % (r["gene"], r["mean_logFC"]) for r in lst[:6])))
print("  %-16s %3d（被多面板共享，不计入归属）" % ("共享标志物", len(shared)))
print("  %-16s %3d（不在任何面板内）" % ("未归属", len(unassigned)))

# ---- 与 WGCNA / 富集结果交集 ----
print("\n" + "-" * 112)
print("与你自己已有产出的交集")
print("-" * 112)
core_set = set(r["gene"] for r in core)
for label, path in [
    ("WGCNA turquoise 模块", os.path.join(D2, "02.WGCNA/WGCNA_top_module_genes.csv")),
    ("WGCNA 全模块", os.path.join(D2, "02.WGCNA/WGCNA_module_genes.csv")),
    ("富集 Hub 基因", os.path.join(D2, "03.enrich/gene_Tumor_vs_Polyp.csv")),
]:
    if not os.path.exists(path):
        print("  %-24s 文件不存在: %s" % (label, os.path.basename(path)))
        continue
    genes = set()
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            for k, v in r.items():
                if v:
                    genes.add(str(v).strip().strip('"'))
    inter = core_set & genes
    print("  %-24s 共 %5d 基因 | 与 Core 共识交集 %4d (%.1f%%)"
          % (label, len(genes), len(inter), 100.0 * len(inter) / max(len(core_set), 1)))

# ---- 落盘 ----
outdir = BASE
with open(os.path.join(outdir, "consensus_all.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene", "n_cohort", "mean_logFC", "median_logFC", "agree", "dir_ratio",
                "z_stouffer", "p_stouffer", "q_stouffer", "p_fisher", "min_p", "cohort_logFC"])
    for r in sorted(rows, key=lambda r: -abs(r["z_stouffer"])):
        w.writerow([r["gene"], r["n_cohort"], round(r["mean_logFC"], 4), round(r["median_logFC"], 4),
                    r["agree"], round(r["dir_ratio"], 2), round(r["z_stouffer"], 3),
                    "%.3e" % r["p_stouffer"], "%.3e" % r["q_stouffer"], "%.3e" % r["p_fisher"],
                    "%.3e" % r["min_p"], "|".join(str(x) for x in r["fcs"])])

with open(os.path.join(outdir, "consensus_core.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["gene", "mean_logFC", "median_logFC", "z_stouffer", "q_stouffer", "cohort_logFC", "change"])
    for r in core:
        w.writerow([r["gene"], round(r["mean_logFC"], 4), round(r["median_logFC"], 4),
                    round(r["z_stouffer"], 3), "%.3e" % r["q_stouffer"],
                    "|".join(str(x) for x in r["fcs"]),
                    "UP" if r["mean_logFC"] > 0 else "DOWN"])

summary = {
    "cohorts": names,
    "n_genes_mergeable": len(rows),
    "core_n": len(core), "core_up": up_core, "core_down": len(core) - up_core,
    "high_n": len(high),
    "core_genes": [r["gene"] for r in core],
    "attribution": {k: [r["gene"] for r in v] for k, v in excl.items()},
    "shared": [[g, o] for g, o in shared],
    "pct_core_in_stroma_activated":
        100.0 * len([g for g in excl.get("间质·激活/CAF", [])]) / max(len(core), 1),
}
json.dump(summary, open(os.path.join(outdir, "consensus_summary.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2, default=str)
print("\n[saved] consensus_all.csv / consensus_core.csv / consensus_summary.json")
print("DONE")
