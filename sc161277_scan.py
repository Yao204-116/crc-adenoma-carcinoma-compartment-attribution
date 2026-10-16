# -*- coding: utf-8 -*-
"""GSE161277 单细胞核查：腺瘤期「分化程序升高」是组成（更多分化细胞）还是状态（每细胞更强）？

设计（与 LCM 分析逻辑对齐，桥接两层数据）：
  组织层信号  T = 全部 QC 细胞平均分化得分        ≈ LCM 测到的东西
  上皮层信号  U = 上皮细胞内平均分化得分           = 单细胞才有的"状态"读数
  上皮占比    P = 上皮细胞数 / 全部 QC 细胞        = 单细胞才有的"组成"读数
  当非上皮细胞该得分为 0 时 T = P × U 精确成立：
      Δlog T = Δlog P（组成项） + Δlog U（状态项）

QC：total UMI ≥ 500 且 pct_mito ≤ 20%
上皮判定：EPCAM/KRT8/18/19/20/CDH1 合计 UMI 占比 ≥ 0.5%（blood 为阴性对照校验阈值）
用法：python sc161277_scan.py
输出：sc161277_scan.json
"""
import gzip, json, os

BASE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(BASE, "sc161277")

SAMPLES = {
    "GSM4904234": ("Patient0", "carcinoma"),
    "GSM4904235": ("Patient1", "adenoma"),
    "GSM4904236": ("Patient1", "carcinoma"),
    "GSM4904237": ("Patient1", "normal"),
    "GSM4904238": ("Patient2", "adenoma"),
    "GSM4904239": ("Patient2", "carcinoma"),
    "GSM4904240": ("Patient2", "normal"),
    "GSM4904241": ("Patient2", "para-cancer"),
    "GSM4904242": ("Patient3", "adenoma"),
    "GSM4904243": ("Patient3", "adenoma"),
    "GSM4904244": ("Patient3", "blood"),
    "GSM4904245": ("Patient3", "carcinoma"),
    "GSM4904246": ("Patient3", "normal"),
}

EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDH1"]
SIGS = {
    "stem":     ["LGR5", "OLFM4", "ASCL2", "SMOC2", "RNF43", "EPHB2", "PROM1"],
    "prolif":   ["MKI67", "TOP2A", "PCNA"],
    "diff_abs": ["ALPI", "FABP1", "FABP2", "CA1", "CA2", "GUCA2A", "GUCA2B", "CDX2", "CEACAM7"],
    "goblet":   ["MUC2", "TFF3", "FCGBP", "ZG16", "SPDEF", "ATOH1"],
    "paneth":   ["DEFA5", "DEFA6", "LYZ", "REG3G", "PLA2G2A", "ITLN1"],
}
IMMUNE = ["PTPRC", "CD3D", "CD8A", "MS4A1", "LYZ", "LST1", "NKG7"]
STROMA = ["COL1A1", "DCN", "VIM", "ACTA2", "TAGLN"]

WANT = sorted(set(EPI + [g for v in SIGS.values() for g in v] + IMMUNE + STROMA))
EPI_MIN_FRAC = 0.005
UMI_MIN = 500
MT_MAX = 0.20


def read_features(path):
    """返回 (row→symbol dict, symbol→[rows])"""
    r2s, s2r = {}, {}
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            p = line.rstrip("\n").split("\t")
            g = p[1] if len(p) >= 2 else (p[0] if p else "")
            r2s[i] = g
            s2r.setdefault(g, []).append(i)
    return r2s, s2r


def scan_sample(gsm, stem):
    feat_path = os.path.join(D, stem + "_features.tsv.gz")
    mtx_path = os.path.join(D, stem + "_matrix.mtx.gz")
    r2s, s2r = read_features(feat_path)
    ncell = 0
    with gzip.open(os.path.join(D, stem + "_barcodes.tsv.gz"), "rt") as f:
        for _ in f:
            ncell += 1

    want_rows = {}
    for g in WANT:
        for r in s2r.get(g, []):
            want_rows[r] = g
    mt_rows = set(r for r, g in r2s.items() if g.startswith("MT-"))
    epi_rows = set(r for g in EPI for r in s2r.get(g, []))

    cell_total = [0.0] * ncell
    cell_mt = [0.0] * ncell
    epi_umis = [0.0] * ncell
    gcount = {g: [0.0] * ncell for g in WANT}

    nline = 0
    with gzip.open(mtx_path, "rt") as f:
        header = f.readline()
        while header.startswith("%"):
            header = f.readline()
        nr, nc, nnz = (int(x) for x in header.split())
        assert nc == ncell, "matrix/barcode 不一致 %s: mtx=%d bc=%d" % (gsm, nc, ncell)
        for line in f:
            p = line.split()
            if len(p) != 3:
                continue
            r = int(p[0]) - 1
            c = int(p[1]) - 1
            v = float(p[2])
            nline += 1
            cell_total[c] += v
            if r in mt_rows:
                cell_mt[c] += v
            if r in epi_rows:
                epi_umis[c] += v
            g = want_rows.get(r)
            if g is not None:
                gcount[g][c] += v
    assert nline == nnz, "nnz 不一致 %s: %d vs %d" % (gsm, nline, nnz)
    return barcodes_of(stem), gcount, cell_total, cell_mt, epi_umis, ncell


def barcodes_of(stem):
    with gzip.open(os.path.join(D, stem + "_barcodes.tsv.gz"), "rt") as f:
        return [l.strip() for l in f if l.strip()]


def mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def analyse(stem, gsm, stage):
    bc, gcount, cell_total, cell_mt, epi_umis, ncell = scan_sample(gsm, stem)
    # QC + 上皮判定 + 归一化得分
    keep = []
    for c in range(ncell):
        t = cell_total[c]
        if t < UMI_MIN:
            continue
        if cell_mt[c] / t > MT_MAX:
            continue
        keep.append(c)
    epi_idx, non_idx = [], []
    for c in keep:
        (epi_idx if epi_umis[c] / cell_total[c] >= EPI_MIN_FRAC else non_idx).append(c)

    # CP10K + log1p
    def norm_expr(g, c):
        return __import__("math").log1p(gcount[g][c] / cell_total[c] * 1e4)

    def sig_score(cells, genes):
        if not cells:
            return 0.0
        return mean([mean([norm_expr(g, c) for g in genes]) for c in cells])

    def frac_pos(cells, genes, thr=0.0):
        """细胞内至少一个基因 raw>0 的比例（检出的宽松版）"""
        if not cells:
            return 0.0
        n = 0
        for c in cells:
            if any(gcount[g][c] > 0 for g in genes):
                n += 1
        return n / len(cells)

    out = {
        "gsm": gsm, "stage": stage,
        "n_cells_raw": ncell, "n_cells_qc": len(keep),
        "pct_epi_of_qc": len(epi_idx) / len(keep) if keep else 0.0,
        "median_umi": sorted(cell_total[c] for c in keep)[len(keep) // 2] if keep else 0,
    }
    # 三层读数：对每个签名
    for name, genes in SIGS.items():
        out[name] = {
            "T_all_cells": sig_score(keep, genes),          # 组织层 ≈ LCM
            "U_epi_cells": sig_score(epi_idx, genes),        # 上皮内状态
            "T_non_epi": sig_score(non_idx, genes),          # 应接近免疫/基质信号
            "frac_epi_pos": frac_pos(epi_idx, genes),        # 上皮内检出比例（组成 proxy）
        }
    # 免疫/基质校验
    out["immune_U_epi"] = sig_score(epi_idx, IMMUNE)
    out["immune_T_all"] = sig_score(keep, IMMUNE)
    out["stroma_U_epi"] = sig_score(epi_idx, STROMA)
    out["stroma_T_all"] = sig_score(keep, STROMA)
    out["epi_core_U_epi"] = sig_score(epi_idx, EPI)
    out["epi_core_T_all"] = sig_score(keep, EPI)
    # 潘氏基因逐一（ad 峰是 DEFA6/HD5 驱动的，要看单基因）
    out["genes_U_epi"] = {g: sig_score(epi_idx, [g]) for g in
                          ["DEFA5", "DEFA6", "LYZ", "REG3G", "PLA2G2A", "ITLN1",
                           "MUC2", "TFF3", "ALPI", "FABP1", "CA2", "GUCA2A",
                           "LGR5", "OLFM4", "MKI67", "EPCAM", "KRT20", "CDX2"]}
    # per-patient/stage 汇总用
    out["_epi_idx_n"] = len(epi_idx)
    return out


def main():
    results = []
    for gsm, (pat, stage) in sorted(SAMPLES.items()):
        stem = None
        for fn in os.listdir(D):
            if fn.startswith(gsm) and fn.endswith("_matrix.mtx.gz"):
                stem = fn[: -len("_matrix.mtx.gz")]
        if stem is None or os.path.getsize(os.path.join(D, stem + "_matrix.mtx.gz")) < 1_000_000:
            print("跳过（matrix 未下载/不完整）:", gsm, flush=True)
            continue
        print("处理", gsm, pat, stage, flush=True)
        try:
            r = analyse(stem, gsm, stage)
        except (gzip.BadGzipFile, EOFError, AssertionError) as e:
            print("  !! 跳过损坏/不完整矩阵 %s: %r" % (gsm, e), flush=True)
            continue
        r["patient"] = pat
        results.append(r)
    json.dump(results, open(os.path.join(BASE, "sc161277_scan.json"), "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)

    # 摘要打印
    print("\n%-8s %-11s %7s %7s %7s | %7s %7s %7s" %
          ("patient", "stage", "raw", "qc", "epi%", "T_stem", "U_stem", "T_prolif"))
    for r in results:
        print("%-8s %-11s %7d %7d %6.1f%% | %7.3f %7.3f %7.3f" % (
            r["patient"], r["stage"], r["n_cells_raw"], r["n_cells_qc"],
            r["pct_epi_of_qc"] * 100,
            r["stem"]["T_all_cells"], r["stem"]["U_epi_cells"], r["prolif"]["T_all_cells"]))


if __name__ == "__main__":
    main()
