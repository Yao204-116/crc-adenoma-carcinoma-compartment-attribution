# -*- coding: utf-8 -*-
"""
诊断 GSE15960（LCM 纯上皮队列）里关键基因的分组均值，用于交叉验证 sig_validate 的结论。
动机：sig_validate 报出「下调分量在 LCM 里 腺瘤 > 正常 > 癌」（非单调），这在生物学上需要核实，
     不能直接写进报告。这里把原始 log2 均值按组打出来，做人工判读。
同时输出 AUC 的方向说明（AUC<0.5 = 反向完美区分，不是"无区分力"）。
"""
import csv
import gzip
import math
import os
import statistics
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "cohort_hunt", "geo_cache")
UA = {"User-Agent": "Mozilla/5.0"}


def load(gse, gpl):
    mf = os.path.join(CACHE, gse + "_series_matrix.txt.gz")
    af = os.path.join(CACHE, gpl + ".annot.gz")
    txt = gzip.open(mf, "rt", encoding="utf-8", errors="replace").read()
    lines = txt.split("\n")

    titles, chars, tbl = [], [], None
    for i, ln in enumerate(lines):
        if ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in ln.split("\t")[1:]])
        elif ln.startswith("!series_matrix_table_begin"):
            tbl = i + 1
            break
    ns = len(titles)
    cbs = [""] * ns
    for row in chars:
        for j in range(min(ns, len(row))):
            cbs[j] += " " + row[j]

    import sys
    sys.path.insert(0, BASE)
    from cohort_hunt2 import classify
    labs = [classify(titles[i], "", cbs[i])[0] for i in range(ns)]

    ann = {}
    atxt = gzip.open(af, "rt", encoding="utf-8", errors="replace").read().split("\n")
    started = False
    for ln in atxt:
        if ln.startswith("!platform_table_begin"):
            started = True
            continue
        if ln.startswith("!platform_table_end"):
            break
        if not started:
            continue
        p = ln.split("\t")
        if len(p) < 3:
            continue
        sym = p[2].strip().strip('"')
        if sym and sym not in ("", "---"):
            ann[p[0].strip().strip('"')] = sym.split("///")[0].strip()

    expr = {}
    rdr = csv.reader(lines[tbl:], delimiter="\t")
    next(rdr, None)
    for row in rdr:
        if not row or row[0].startswith("!"):
            break
        g = ann.get(row[0].strip().strip('"'))
        if not g:
            continue
        vals = []
        for v in row[1:]:
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(None)
        if len(vals) < ns:
            continue
        expr.setdefault(g, []).append(vals)

    # 多探针 → 取跨样本均值最大的探针（与 sig_validate 一致）
    out = {}
    for g, rows in expr.items():
        best = max(rows, key=lambda v: statistics.mean([x for x in v if x is not None]))
        out[g] = best
    # log2
    med = statistics.median([x for v in out.values() for x in v if x is not None])
    if med > 30:
        out = {g: [math.log2(x + 1) if x is not None else None for x in v] for g, v in out.items()}
    return labs, out


def mean_of(vals, labs, want):
    xs = [vals[i] for i, l in enumerate(labs) if l == want and i < len(vals) and vals[i] is not None]
    return statistics.mean(xs) if xs else float("nan")


def main():
    labs, expr = load("GSE15960", "GPL570")
    from collections import Counter
    print("分组构成:", dict(Counter(labs)))
    groups = ["normal", "adenoma", "carcinoma"]
    print("\n%-10s %10s %10s %10s   %s" % ("基因", "正常", "腺瘤", "癌", "腺瘤→癌 变化"))
    print("-" * 78)
    KEY = [
        # 构建期「下调分量」的代表（上皮分化/功能）
        "CLCA1", "FCGBP", "SPINK4", "MUC2", "REG4", "ITLN1", "PIGR", "SI", "GUCA2A", "ZG16",
        # 构建期「上调分量」的代表（间质激活 / 髓系）
        "SPP1", "FAP", "COL1A1", "COL11A1", "THBS2", "ASPN", "GREM1", "SULF1", "CXCL8", "S100A9",
        # 上皮结构基线（用于判断"是否细胞数量变化"）
        "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1",
        # 淋巴细胞基线
        "CD3D", "CD8A", "CD79A", "MS4A1",
    ]
    for g in KEY:
        if g not in expr:
            print("%-10s  (未比对到)" % g)
            continue
        v = expr[g]
        a, b, c = (mean_of(v, labs, k) for k in groups)
        arrow = "↑上升" if c > a else ("↓下降" if c < a else "=")
        print("%-10s %10.2f %10.2f %10.2f   %s (%+.2f)" % (g, a, b, c, arrow, c - a))

    print("\n说明：这是 LCM（激光显微切割）纯上皮样本 → 不含间质/免疫细胞。")
    print("      若某基因在此队列「腺瘤→癌」仍显著变化，说明是上皮内在变化；")
    print("      若几乎不动或反向，说明它在全组织队列里的变化来自间质/免疫细胞的占比（组成效应）。")


if __name__ == "__main__":
    main()
