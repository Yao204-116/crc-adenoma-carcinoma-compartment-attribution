# -*- coding: utf-8 -*-
"""
全局漂移的性质诊断：技术偏移 还是 真实生物学？
================================================================================
问题：GSE15960（LCM 纯上皮）全基因组 Δ(癌-腺瘤) = −0.40，必须判断它的来源，
      因为它决定「是否该校正」，而校正与否会改变 2×2 分解的结论。

三个诊断：
  A. 逐样本全局水平 glev 按组列出 —— 看是「组间整齐分离」还是「组内散乱」
  B. 用「低方差基因」重算 glev（稳定基因是纯技术偏移的最好探针）
     → 若低方差基因给出的组间偏移与全基因组一致，说明偏移是全局性的
  C. 三组排序（正常/腺瘤/癌）—— 单调性 vs 跳变
原理：
  低方差基因（在 18 个样本间几乎不变）几乎不可能有真实的组间生物学差异。
  因此它们之间的组间系统性差异，只能来自技术（RNA 质量、杂交批次、扫描增益）。
用法：python drift_diag.py
"""
import json, math, os, statistics
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
CACHE = os.path.join(HERE, "cohort_hunt", "geo_cache")

_src = open("sig_validate.py", encoding="utf-8").read().split("\n")
_ns = {"__name__": "_svdefs"}
exec("\n".join(_src[:266]), _ns)
parse_series, load_annot, load_entrez_map = _ns["parse_series"], _ns["load_annot"], _ns["load_entrez_map"]

import gzip
T2K = {"tumor": "carcinoma", "adenoma": "adenoma", "ssa": "ssa", "adjacent mucosa": "normal"}


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
    return [T2K.get(t, "other") for t in tis] if tis else None


def load_cohort(gse):
    d = parse_series(gse)
    if d is None:
        return None
    ann = load_annot(d["gpl"]) or load_entrez_map(d["expr"].keys())
    if not ann:
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
    if med and statistics.median(med) > 30:
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
    # 保留原始（log2 后）表达矩阵，用于低方差基因探针
    raw = {g: gene_rows[g][1][:NG] for g in gene_rows}
    return {"zmap": zmap, "raw": raw, "labels": lab, "gpl": d["gpl"]}


COHORTS = ["GSE4183", "GSE117606", "GSE15960"]
KIND = {"GSE4183": "全组织 Affy", "GSE117606": "全组织 Agilent", "GSE15960": "LCM 纯上皮"}
C = {}
for g in COHORTS:
    print("载入 %s ..." % g)
    C[g] = load_cohort(g)

OUT = {}
for gse in COHORTS:
    c = C[gse]
    if not c:
        continue
    lab = c["labels"]
    zmap, raw = c["zmap"], c["raw"]
    NS = len(lab)
    gi = {k: [i for i, l in enumerate(lab) if l == k] for k in ("normal", "adenoma", "carcinoma")}
    print("\n" + "=" * 96)
    print("%s（%s）  正常 %d / 腺瘤 %d / 癌 %d" % (gse, KIND[gse], len(gi["normal"]),
                                                len(gi["adenoma"]), len(gi["carcinoma"])))
    print("=" * 96)

    # ---- A. 逐样本全局水平（全基因组 z 均值）----
    glev = []
    for i in range(NS):
        v = [zmap[g][i] for g in zmap if i < len(zmap[g])]
        glev.append(statistics.mean(v) if v else 0.0)
    print("\n  A. 逐样本全局水平（全基因组 z 均值，样本顺序 = series matrix 列顺序）")
    for k in ("normal", "adenoma", "carcinoma"):
        if not gi[k]:
            continue
        vals = [glev[i] for i in gi[k]]
        print("     %-11s n=%2d  glev: %s   → 组均值 %+0.4f"
              % (k, len(vals), " ".join("%+.3f" % x for x in vals), statistics.mean(vals)))

    # ---- B. 低方差基因重算（稳定基因 = 技术偏移探针）----
    # 只用「表达上有意义（中位 > 全基因中位）」的基因，避免探针噪声干扰
    allmed = [statistics.median([v for v in raw[g] if v is not None]) for g in raw]
    thr = statistics.median(allmed)
    cand = [g for g in raw if statistics.median([v for v in raw[g] if v is not None]) > thr]
    sd_of = {g: statistics.pstdev([v for v in raw[g][:NS] if v is not None]) for g in cand}
    stable = sorted(cand, key=lambda g: sd_of[g])[:3000]
    print("\n  B. 低方差探针（表达高于中位、样本间 SD 最小的 3000 个基因）")
    print("     ↑ 这些基因几乎不可能有真实组间生物学差异 → 它们的组间差异＝技术偏移")
    for k in ("normal", "adenoma", "carcinoma"):
        if not gi[k]:
            continue
        vals = [statistics.mean([raw[g][i] for g in stable if i < len(raw[g]) and raw[g][i] is not None])
                for i in gi[k]]
        print("     %-11s 原始 log2 均值 %8.4f   （逐样本 %s）"
              % (k, statistics.mean(vals), " ".join("%.2f" % x for x in vals)))
    if gi["adenoma"] and gi["carcinoma"]:
        s_ad = statistics.mean([statistics.mean([raw[g][i] for g in stable
                                                if i < len(raw[g]) and raw[g][i] is not None])
                               for i in gi["adenoma"]])
        s_ca = statistics.mean([statistics.mean([raw[g][i] for g in stable
                                                if i < len(raw[g]) and raw[g][i] is not None])
                               for i in gi["carcinoma"]])
        print("     → 稳定基因只差 %+0.4f log2（癌 − 腺瘤）" % (s_ca - s_ad))
    # 低方差基因的 z 全局水平
    zstable = {g: zmap[g] for g in stable if g in zmap}
    gl2 = []
    for i in range(NS):
        v = [zstable[g][i] for g in zstable]
        gl2.append(statistics.mean(v) if v else 0.0)
    print("     稳定基因的 z 全局水平：", " ".join("%+.3f" % x for x in gl2))
    for k in ("normal", "adenoma", "carcinoma"):
        if gi[k]:
            print("     %-11s → %+0.4f" % (k, statistics.mean([gl2[i] for i in gi[k]])))

    # ---- C. 单调性 ----
    print("\n  C. 三组排序（判断是「单调梯度」还是「单组跳变」）")
    order = [("normal", "正常"), ("adenoma", "腺瘤"), ("carcinoma", "癌")]
    gm = {k: statistics.mean([glev[i] for i in gi[k]]) for k, _ in order if gi[k]}
    print("     glev:  " + "  →  ".join("%s %+0.3f" % (nm, gm[k]) for k, nm in order if k in gm))
    mm = {}
    for k, _ in order:
        if gi[k]:
            mm[k] = 100.0 * sum(1 for g in zmap
                                if statistics.mean([zmap[g][i] for i in gi[k] if i < len(zmap[g])]) > 0) / len(zmap)
    print("     上调基因占比:  " + "  →  ".join("%s %.1f%%" % (nm, mm[k]) for k, nm in order if k in mm))
    OUT[gse] = {"glev": glev, "stable_glev_mean": {k: statistics.mean([gl2[i] for i in gi[k]])
                                                   for k, _ in order if gi[k]},
                "group_glev": gm, "pct_up": mm}

json.dump(OUT, open("drift_diag.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\n→ drift_diag.json")
print("DONE")
