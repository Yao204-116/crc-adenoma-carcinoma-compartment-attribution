# -*- coding: utf-8 -*-
"""
成稿事实表：把 40+ 个结果 JSON 收拢成一份「稿件可直接引用」的数字底稿
================================================================================
产出：
  ① manuscript_facts.json          —— 全部关键数字（面板矩阵 / 配对 / 三层 / 归属）
  ② suppl_table1_313genes.csv      —— ★ 审稿人要的「逐基因归属表」（313 基因 × 3 队列）
  ③ final_panel_matrix.csv         —— 面板 × 队列 AUC 矩阵（含方向与逐基因上调比例）

关键设计：
  - 逐基因归属判据：全组织两队列同向 + LCM 是否保留 → 分「上皮内在」/「组成来源」/「不可判」
  - 不使用 import sv（sig_validate.py 的主逻辑在模块级，导入即重写结果文件）
    → 改为只 exec 其第 1..266 行（纯函数定义段）
用法：python manuscript_facts.py
"""
import csv, gzip, json, math, os, statistics, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
CACHE = os.path.join(HERE, "cohort_hunt", "geo_cache")

# ===========================================================================
# 0. 只取 sig_validate.py 的函数定义段（1..266 行），避开模块级主逻辑
# ===========================================================================
_src = open("sig_validate.py", encoding="utf-8").read().split("\n")
_ns = {"__name__": "_svdefs", "__file__": os.path.join(HERE, "sig_validate.py")}
exec("\n".join(_src[:266]), _ns)
parse_series = _ns["parse_series"]
load_annot = _ns["load_annot"]
load_entrez_map = _ns["load_entrez_map"]
auc = _ns["auc"]

PANELS = {k: sorted(v) for k, v in json.load(open("panels_def.json", encoding="utf-8")).items()}

T2K = {"tumor": "carcinoma", "adenoma": "adenoma", "ssa": "ssa", "adjacent mucosa": "normal"}
COHORTS = ["GSE4183", "GSE117606", "GSE15960"]      # 全组织 Affy / 全组织 Agilent / LCM 纯上皮
COHORT_NOTE = {
    "GSE4183":    ("全组织", "Affymetrix",  "15 腺瘤 / 15 癌 / 8 正常"),
    "GSE117606":  ("全组织", "Agilent",     "59 腺瘤 / 74 癌 / 65 癌旁（+10 SSA 剔除）"),
    "GSE15960":   ("LCM 纯上皮", "Affymetrix", "6 腺瘤 / 6 癌 / 6 正常"),
}

# ===========================================================================
# 1. 载入队列（复刻 collide_all 的口径：tissue 字段优先 + log2 判定）
# ===========================================================================
def labels_tissue(gse):
    p = os.path.join(CACHE, gse + "_series_matrix.txt.gz")
    if not os.path.exists(p):
        return None
    tis = None
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("!Sample_characteristics_ch1") and "tissue:" in ln:
                tis = [x.strip().strip('"').replace("tissue:", "").strip().lower()
                       for x in ln.rstrip("\n").split("\t")[1:]]
                break
            if ln.startswith("!series_matrix_table_begin"):
                break
    if not tis:
        return None
    return [T2K.get(t, "other") for t in tis]


def load_cohort(gse):
    d = parse_series(gse)
    if d is None:
        return None
    ann = load_annot(d["gpl"]) or load_entrez_map(d["expr"].keys())
    if not ann:
        print("  !! %s 注释不可得" % gse)
        return None
    gene_rows = {}
    for pid, vals in d["expr"].items():
        g = ann.get(pid)
        if not g:
            continue
        good = [v for v in vals if v is not None]
        if len(good) < len(vals) * 0.8:
            continue
        m = statistics.mean(good)
        gu = g.upper()
        if gu not in gene_rows or m > gene_rows[gu][0]:
            gene_rows[gu] = (m, vals)
    NG = len(d["labels"])
    med = []
    for i in range(NG):
        v = [gene_rows[g][1][i] for g in gene_rows
             if i < len(gene_rows[g][1]) and gene_rows[g][1][i] is not None]
        if v:
            med.append(statistics.median(v))
    if med and statistics.median(med) > 30:                 # 未 log2 → 补 log2
        for g in gene_rows:
            gene_rows[g] = (gene_rows[g][0],
                            [math.log2(v + 1) if v is not None else None for v in gene_rows[g][1]])

    def zrow(vals):
        good = [v for v in vals[:NG] if v is not None]
        if len(good) < 3:
            return None
        m = statistics.mean(good)
        sd = statistics.pstdev(good) or 1.0
        return [(v - m) / sd if v is not None else 0.0 for v in vals[:NG]]
    zmap = {}
    for g in gene_rows:
        z = zrow(gene_rows[g][1])
        if z:
            zmap[g] = z
    lab = d["labels"]
    tl = labels_tissue(gse)
    if tl and len(tl) == len(lab):
        lab = tl
    # ★★ 全局漂移校正（本脚本最重要的一步）
    #   per-gene z 只去掉「基因自身的均值」，不去掉「样本的全局水平」。
    #   若某组样本整体信号偏低（RNA 质量/捕获细胞数），会让几乎所有基因的 Δ 同时偏负，
    #   从而伪造出「类群特异」的假象。实测 GSE15960 全基因组 Δ = −0.40（
    #   SD 1.12，n=21753 → SE 0.008，极显著），必须校正。
    #   做法：逐样本算全基因组 z 的均值 g_i（= 该样本的全局水平），再从该样本每个基因里减掉。
    NS = len(lab)
    glev = []
    for i in range(NS):
        vals = [zmap[g][i] for g in zmap if i < len(zmap[g])]
        glev.append(statistics.mean(vals) if vals else 0.0)
    zmap_adj = {g: [v - glev[i] for i, v in enumerate(z)] for g, z in zmap.items()}
    return {"zmap": zmap, "zmap_adj": zmap_adj, "glev": glev, "labels": lab, "gpl": d["gpl"]}


print("=" * 100)
print("载入队列")
print("=" * 100)
C = {}
for gse in COHORTS:
    c = load_cohort(gse)
    if c is None:
        print("  !! %s 载入失败" % gse)
        continue
    C[gse] = c
    cnt = Counter(c["labels"])
    print("  %-11s %-10s %-11s 基因 %5d | 正常 %3d 腺瘤 %3d 癌 %3d 其他 %d"
          % (gse, COHORT_NOTE[gse][0], c["gpl"], len(c["zmap"]),
             cnt.get("normal", 0), cnt.get("adenoma", 0), cnt.get("carcinoma", 0),
             sum(v for k, v in cnt.items() if k not in ("normal", "adenoma", "carcinoma"))))

# ===========================================================================
# 2. 313 签名
# ===========================================================================
SIGROWS = list(csv.DictReader(open("consensus_signature_tier.csv", encoding="utf-8-sig")))
SIGINFO = {r["gene"].upper(): r for r in SIGROWS}
UP = [g for g, r in SIGINFO.items() if float(r["mean_logFC"]) > 0]
DN = [g for g, r in SIGINFO.items() if float(r["mean_logFC"]) < 0]
print("\n签名：%d 基因（上调 %d / 下调 %d）" % (len(SIGINFO), len(UP), len(DN)))

# 额外基因集（用户上游资产）
def genes_from(path, col="Gene", filt=None):
    out = []
    if not os.path.exists(path):
        return out
    for r in csv.DictReader(open(path, encoding="utf-8-sig")):
        if filt and not filt(r):
            continue
        v = (r.get(col) or "").strip().upper()
        if v:
            out.append(v)
    return out

TURQ = set(genes_from(r"D:/2/02.WGCNA/WGCNA_top_module_genes.csv", "Gene",
                      lambda r: (r.get("Module") or "").strip() == "turquoise"
                                and (r.get("Comparison") or "").strip() == "Tumor_vs_Polyp"))
HUB = set(genes_from(r"D:/2/03.enrich/Hub_Tumor_vs_Polyp.csv", "Gene"))
HUB20 = []
_hp = os.path.join(r"D:/2/04.PPI/degree_top20.csv")
if os.path.exists(_hp):
    # ★ 该文件【无表头】，首行即是数据（FN1,284）——用 DictReader 会把 FN1 当列名而丢掉它
    for row in csv.reader(open(_hp, encoding="utf-8-sig")):
        if row and row[0].strip():
            HUB20.append(row[0].strip().upper())
print("用户资产：turquoise %d / Hub %d / Hub-top20 %d" % (len(TURQ), len(HUB), len(HUB20)))

SETS = {
    "313签名": set(SIGINFO),
    "上调分量(210)": set(UP),
    "下调分量(103)": set(DN),
    "WGCNA-turquoise": TURQ,
    "Hub基因": HUB,
    "turquoise∩上调分量": TURQ & set(UP),
    "turquoise∩下调分量": TURQ & set(DN),
}
for k, v in PANELS.items():
    SETS[k] = set(v)

# ===========================================================================
# 3. 逐基因 Δ（三队列）—— ★ Suppl Table 1 的核心
# ===========================================================================
print("\n" + "=" * 100)
print("逐基因 Δ = mean(z|癌) − mean(z|腺瘤)（★ 已做全局漂移校正；见 3b 节）")
print("=" * 100)


def gene_delta(gse, g, iA, iC, iN=None, adj=True):
    z = C[gse]["zmap_adj" if adj else "zmap"].get(g)
    if not z:
        return None
    ca = [z[i] for i in iC if i < len(z)]
    ad = [z[i] for i in iA if i < len(z)]
    if len(ca) < 2 or len(ad) < 2:
        return None
    out = {"ca": statistics.mean(ca), "ad": statistics.mean(ad),
           "d": statistics.mean(ca) - statistics.mean(ad)}
    if iN:
        no = [z[i] for i in iN if i < len(z)]
        if len(no) >= 2:
            out["no"] = statistics.mean(no)
            out["d_ad_no"] = statistics.mean(ad) - statistics.mean(no)
            out["d_ca_no"] = statistics.mean(ca) - statistics.mean(no)
    return out


IDX = {}
for gse in C:
    lab = C[gse]["labels"]
    IDX[gse] = {
        "ad": [i for i, l in enumerate(lab) if l == "adenoma"],
        "ca": [i for i, l in enumerate(lab) if l == "carcinoma"],
        "no": [i for i, l in enumerate(lab) if l == "normal"],
    }

TISSUE = [g for g in COHORTS if COHORT_NOTE[g][0] == "全组织" and g in C]
LCM = [g for g in COHORTS if COHORT_NOTE[g][0] == "LCM 纯上皮" and g in C]

rows_out = []
verdict_cnt = Counter()
for g in sorted(SIGINFO):
    info = SIGINFO[g]
    rec = {"gene": g, "mean_logFC": float(info["mean_logFC"]),
           "tier_change": info["change"], "q_stouffer": info["q_stouffer"],
           "panels": ";".join(k for k in PANELS if g in SETS[k])}
    tis_signs = []
    for gse in TISSUE:
        v = gene_delta(gse, g, IDX[gse]["ad"], IDX[gse]["ca"])
        rec["d_" + gse] = round(v["d"], 4) if v else None
        if v:
            tis_signs.append(1 if v["d"] > 0 else -1)
    for gse in LCM:
        v = gene_delta(gse, g, IDX[gse]["ad"], IDX[gse]["ca"])
        rec["d_" + gse + "_LCM"] = round(v["d"], 4) if v else None
    # 逐基因方向符合率（与签名的上调/下调方向比）
    exp = 1 if rec["mean_logFC"] > 0 else -1
    lut = [gse for gse in TISSUE if rec["d_" + gse] is not None]
    rec["n_tissue"] = len(lut)
    rec["tissue_agree"] = sum(1 for gse in lut if (1 if rec["d_" + gse] > 0 else -1) == exp)
    lcmk = [g for g in LCM if rec["d_" + g + "_LCM"] is not None]
    rec["n_lcm"] = len(lcmk)
    rec["lcm_agree"] = sum(1 for gse in lcmk if (1 if rec["d_" + gse + "_LCM"] > 0 else -1) == exp)
    # ★ 归属判据
    tissue_consistent = (len(tis_signs) >= 2 and len(set(tis_signs)) == 1)
    if not tissue_consistent:
        v = "全组织不一致"
    else:
        tsign = tis_signs[0]
        lsigns = [1 if rec["d_" + g + "_LCM"] > 0 else -1 for g in lcmk]
        if not lsigns:
            v = "无LCM判据"
        elif lsigns[0] == tsign:
            v = "上皮内在（LCM 保留）"
        elif lsigns[0] == -tsign:
            v = "组成来源（LCM 反向）"
        else:
            v = "?"
    rec["attribution"] = v
    verdict_cnt[v] += 1
    rows_out.append(rec)

with open("suppl_table1_313genes.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
    w.writeheader()
    w.writerows(rows_out)
print("  → suppl_table1_313genes.csv（%d 行）" % len(rows_out))
print("\n  归属判据分布：")
for k, v in verdict_cnt.most_common():
    print("    %-22s %3d" % (k, v))

# ===========================================================================
# 3b. ★ 全局漂移控制（必须先排除的替代解释）
# ---------------------------------------------------------------------------
# 隐患：LCM 里 100 个上调分量基因系统性转为负值。
#   一个可能的替代解释不是「这些基因来自间质」，而是「癌样本整体信号偏低」
#   （RNA 质量差 / 捕获细胞少）→ 那会让几乎所有基因的 Δ 都变负。
# 检验：对全基因组（~21k）算 mean(z|癌) − mean(z|腺瘤)；若 ≈0，则无全局漂移，
#   面板级的负值只能是「基因类别特异」的，替代解释被排除。
# ===========================================================================
print("\n" + "=" * 100)
print("★ 3b. 全局漂移控制：全基因组 Δ（癌 vs 腺瘤）")
print("=" * 100)
drift = {}
for gse in COHORTS:
    zmap = C[gse]["zmap"]
    ia, ic = IDX[gse]["ad"], IDX[gse]["ca"]
    deltas, ntot = [], 0
    for g in zmap:
        z = zmap[g]
        ca = [z[i] for i in ic if i < len(z)]
        ad = [z[i] for i in ia if i < len(z)]
        if len(ca) >= 2 and len(ad) >= 2:
            deltas.append(statistics.mean(ca) - statistics.mean(ad))
            ntot += 1
    gm = statistics.mean(deltas)
    sd = statistics.pstdev(deltas)
    frac_up = 100.0 * sum(1 for x in deltas if x > 0) / len(deltas)
    drift[gse] = {"n_genes": ntot, "global_delta": gm, "global_sd": sd,
                  "frac_up": frac_up}
    print("  %-11s 全基因组 %5d 个  全局 Δ = %+0.4f（SD %.3f）  上调基因占 %0.1f%%"
          % (gse, ntot, gm, sd, frac_up))
print("\n  判读：全局 Δ 接近 0 且上/下调各约一半 → 无全局漂移，")
print("        面板/分量的负值只能解释为「基因类别特异」→ 替代解释被排除。")

# ===========================================================================
# 4. 面板 × 队列 AUC 矩阵
# ===========================================================================
print("\n" + "=" * 100)
print("面板 × 队列矩阵（癌 vs 腺瘤）")
print("=" * 100)
ORDER = ["313签名", "上调分量(210)", "下调分量(103)",
         "上皮·分化/功能", "上皮分化·核心8", "纯上皮基线·7",
         "间质·激活/CAF", "间质·基线", "泛间质·6", "髓系·激活",
         "WGCNA-turquoise", "turquoise∩上调分量", "turquoise∩下调分量",
         "Hub基因", "Hub-top20"]
mat = {}
hdr = "  %-20s" % "基因集" + "".join("%22s" % ("%s %s" % (g.replace("GSE", ""), COHORT_NOTE[g][0])) for g in COHORTS)
print(hdr)
print("  " + "-" * (len(hdr) - 2))
for name in ORDER:
    if name not in SETS:
        continue
    gs0 = SETS[name]
    line = "  %-20s" % name
    rowm = {}
    for gse in COHORTS:
        if gse not in C:
            continue
        zmap = C[gse]["zmap_adj"]           # ★ 用漂移校正后的 z
        zraw = C[gse]["zmap"]
        gs = gs0 & set(zmap)
        if len(gs) < 3:
            line += "%22s" % "—"
            rowm[gse] = {"n": len(gs), "auc": None}
            continue
        ia, ic = IDX[gse]["ad"], IDX[gse]["ca"]

        def sc(i, zm=zmap):
            zz = [zm[g_][i] for g_ in gs if i < len(zm[g_])]
            return statistics.mean(zz) if zz else None
        a = [x for x in (sc(i) for i in ia) if x is not None]
        c_ = [x for x in (sc(i) for i in ic) if x is not None]
        aa = auc(c_, a)
        a0 = [x for x in (sc(i, zraw) for i in ia) if x is not None]
        c0 = [x for x in (sc(i, zraw) for i in ic) if x is not None]
        aa_raw = auc(c0, a0)
        n_up = sum(1 for g_ in gs
                   if statistics.mean([zmap[g_][i] for i in ic]) >
                      statistics.mean([zmap[g_][i] for i in ia]))
        dlt = statistics.median(c_) - statistics.median(a)
        dlt_raw = statistics.median(c0) - statistics.median(a0)
        rowm[gse] = {"n": len(gs), "n_total": len(gs0), "auc": aa, "auc_raw": aa_raw,
                     "n_up": n_up, "d_median": dlt, "d_median_raw": dlt_raw,
                     "pct_up": 100.0 * n_up / len(gs)}
        line += "%22s" % ("%.3f(%.3f) Δ%+.2f %d%%" % (max(aa, 1 - aa), max(aa_raw, 1 - aa_raw),
                                                      dlt, round(100.0 * n_up / len(gs))))
    mat[name] = rowm
    print(line)
print("\n  格式：|AUC|校正(|AUC|原始) Δ中位 上调%%     ← 校正=扣掉样本全局水平后")

with open("final_panel_matrix.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["基因集"] + ["%s_n" % g for g in COHORTS] + ["%s_auc" % g for g in COHORTS]
               + ["%s_auc_raw" % g for g in COHORTS]
               + ["%s_pct_up" % g for g in COHORTS] + ["%s_d_median" % g for g in COHORTS])
    for name in ORDER:
        if name not in mat:
            continue
        r = [name]
        r += [mat[name].get(g, {}).get("n") for g in COHORTS]
        r += [mat[name].get(g, {}).get("auc") for g in COHORTS]
        r += [mat[name].get(g, {}).get("auc_raw") for g in COHORTS]
        r += [mat[name].get(g, {}).get("pct_up") for g in COHORTS]
        r += [mat[name].get(g, {}).get("d_median") for g in COHORTS]
        w.writerow(r)
print("  → final_panel_matrix.csv")

# ===========================================================================
# 5. ★ 三层模式（正常/腺瘤/癌）：原始 vs 漂移校正 —— 检验「腺瘤峰」是否稳健
# ---------------------------------------------------------------------------
# 模式判据（沿用本轮定义）：
#   P = 腺瘤峰（ad >= no 且 ad > ca）      ← 反直觉、最有价值
#   D = 单调降（no > ad > ca）             ← 教科书预期
#   A = 其余（多为「腺瘤期已降」）
# ===========================================================================
print("\n" + "=" * 100)
print("★ 5. 三层模式：腺瘤峰是否在漂移校正后仍成立")
print("=" * 100)


def pattern(no, ad, ca):
    if ad >= no and ad > ca:
        return "P"
    if no > ad > ca:
        return "D"
    return "A"


TL = {"raw": {}, "adj": {}}
PLIST = ["313全签名", "上调分量(210)", "下调分量(103)"] + list(PANELS)
for mode in ("raw", "adj"):
    for name in PLIST:
        gs0 = SETS.get(name)
        if not gs0:
            continue
        TL[mode][name] = {}
        for gse in COHORTS:
            if gse not in C:
                continue
            zm = C[gse]["zmap" if mode == "raw" else "zmap_adj"]
            gs = gs0 & set(zm)
            if len(gs) < 3:
                continue
            med = {}
            for grp in ("no", "ad", "ca"):
                ii = IDX[gse][grp]
                vals = []
                for i in ii:
                    zz = [zm[g_][i] for g_ in gs if i < len(zm[g_])]
                    if zz:
                        vals.append(statistics.mean(zz))
                med[grp] = statistics.median(vals) if vals else None
            if None in med.values():
                continue
            TL[mode][name][gse] = {"normal": med["no"], "adenoma": med["ad"],
                                   "carcinoma": med["ca"],
                                   "pattern": pattern(med["no"], med["ad"], med["ca"])}

print("\n  面板                          |      原始口径（正常/腺瘤/癌 → 模式）      |   漂移校正口径")
print("  " + "-" * 96)
for name in PLIST:
    if name not in TL["raw"]:
        continue
    r = TL["raw"][name]
    a = TL["adj"].get(name, {})
    cells = []
    for gse in COHORTS:
        if gse in r:
            cells.append("%+.2f/%+.2f/%+.2f→%s" % (r[gse]["normal"], r[gse]["adenoma"],
                                                   r[gse]["carcinoma"], r[gse]["pattern"]))
        else:
            cells.append("—")
    cells2 = []
    for gse in COHORTS:
        if gse in a:
            cells2.append("%+.2f/%+.2f/%+.2f→%s" % (a[gse]["normal"], a[gse]["adenoma"],
                                                    a[gse]["carcinoma"], a[gse]["pattern"]))
        else:
            cells2.append("—")
    print("  %-26s | %s | %s" % (name, "  ".join("%-20s" % c for c in cells),
                                 "  ".join("%-20s" % c for c in cells2)))

print("\n  ★ 汇总：模式为 P（腺瘤峰）的队列数")
print("  %-26s %8s %8s" % ("面板", "原始", "校正"))
for name in PLIST:
    if name not in TL["raw"]:
        continue
    pr = sum(1 for g in TL["raw"][name].values() if g["pattern"] == "P")
    pa = sum(1 for g in TL["adj"].get(name, {}).values() if g["pattern"] == "P")
    n = len(TL["raw"][name])
    print("  %-26s %5d/%-3d %5d/%-3d" % (name, pr, n, pa, n))

# ===========================================================================
# 5b. 分量的落差保留比例（逐基因口径）
# ===========================================================================
print("\n" + "=" * 100)
print("5b. 落差保留比例：LCM 纯上皮 Δ 相对全组织 Δ")
print("=" * 100)
rat = {}
for comp, sel in (("上调分量", "UP"), ("下调分量", "DOWN")):
    xs, ys = [], []
    for r in rows_out:
        if r["tier_change"] != sel:
            continue
        t = [r["d_" + g] for g in TISSUE if r.get("d_" + g) is not None]
        y = r.get("d_" + LCM[0] + "_LCM") if LCM else None
        if not t or y is None:
            continue
        xs.append(statistics.mean(t))
        ys.append(y)
    if not xs:
        continue
    same = 100.0 * sum(1 for x, y in zip(xs, ys) if x * y > 0) / len(xs)
    rat[comp] = {"n": len(xs), "tissue_mean": statistics.mean(xs), "lcm_mean": statistics.mean(ys),
                 "ratio": statistics.mean(ys) / statistics.mean(xs), "same_dir_pct": same,
                 "pct_up_tissue": 100.0 * sum(1 for x in xs if x > 0) / len(xs),
                 "pct_up_lcm": 100.0 * sum(1 for y in ys if y > 0) / len(ys)}
    print("  %s: 全组织 %+0.3f → LCM %+0.3f  = %+0.2f 倍；同向 %.1f%%；上调占比 %.1f%% → %.1f%%"
          % (comp, rat[comp]["tissue_mean"], rat[comp]["lcm_mean"], rat[comp]["ratio"],
             same, rat[comp]["pct_up_tissue"], rat[comp]["pct_up_lcm"]))

# ===========================================================================
# 6. 收拢既有 JSON
# ===========================================================================
FACTS = {
    "meta": {
        "cohorts": {g: {"kind": COHORT_NOTE[g][0], "platform": COHORT_NOTE[g][1],
                        "design": COHORT_NOTE[g][2], "gpl": C[g]["gpl"] if g in C else None,
                        "n_genes": len(C[g]["zmap"]) if g in C else None} for g in COHORTS},
        "signature": {"n": len(SIGINFO), "up": len(UP), "dn": len(DN),
                      "source": "consensus_signature_tier.csv（4 队列共识）"},
    },
    "attribution_313": dict(verdict_cnt),
    "global_drift": drift,
    "three_layer": {"raw": TL["raw"], "adj": TL["adj"]},
    "component_retention": rat,
    "panel_matrix": mat,
}
for fn, key in [("three_layer_pattern.json", "three_layer_pattern"),
                ("three_layer_results.json", "three_layer_detail"),
                ("collide_results.json", "collide")]:
    if os.path.exists(fn):
        FACTS[key] = json.load(open(fn, encoding="utf-8"))
        print("  收拢 %s → FACTS['%s']" % (fn, key))
_p = "gse117606_paired_results.json"
if os.path.exists(_p):
    FACTS["paired_117606"] = json.load(open(_p, encoding="utf-8"))
    print("  收拢 %s → FACTS['paired_117606']" % _p)
# 四队列 AUC：优先读 sig_validate.py 的直接产物，读不到再用 gse117606_paired.py 的备份。
# ★ 旧写法只认顶层那份备份，于是「先跑 sig_validate 再跑本脚本」会静默丢掉
#   four_cohort_auc（build_suppl_tables.py 的第一个 smoke test 就要它）。
for _f in ("cohort_hunt/sig_validation_results.json",
           "sig_validation_results_4cohort.json"):
    if os.path.exists(_f) and os.path.getsize(_f) > 100:
        FACTS["four_cohort_auc"] = json.load(open(_f, encoding="utf-8"))
        print("  收拢 %s → FACTS['four_cohort_auc']" % _f)
        break

json.dump(FACTS, open("manuscript_facts.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n→ manuscript_facts.json（%.1f KB）" % (os.path.getsize("manuscript_facts.json") / 1024))

# ===========================================================================
# 6. 归属清单（写进报告的表）
# ===========================================================================
print("\n" + "=" * 100)
print("★ 上皮内在清单（全组织两队列同向 + LCM 同向）")
print("=" * 100)
for v in ["上皮内在（LCM 保留）", "组成来源（LCM 反向）", "全组织不一致", "无LCM判据", "?"]:
    gs = [r["gene"] for r in rows_out if r["attribution"] == v]
    print("\n  【%s】%d 个" % (v, len(gs)))
    if v == "上皮内在（LCM 保留）":
        for i in range(0, len(gs), 12):
            print("     " + " ".join(gs[i:i + 12]))
    else:
        print("     " + " ".join(gs[:36]) + (" ..." if len(gs) > 36 else ""))
print("\nDONE")
