# -*- coding: utf-8 -*-
"""补充分析 ③：GSE161277 同接型双细胞检测 + 过滤后敏感性分解（B 类补分析）

目的：排除"单细胞腺瘤峰由双细胞（上皮-间质/髓系同接）造成"的指控。

方法（Scrublet/scDblFinder 式，逐样本独立做）：
  1. QC 与 sc161277_scan.py 完全一致（UMI ≥ 500 且线粒体 ≤ 20%），逐样本断言 QC 数一致
  2. 基因过滤（≥3 细胞表达）→ 中位数 CP 归一 + log1p → 高变基因 top 3000 → PCA 50 维
  3. 模拟双细胞 = 随机细胞对的 raw 计数相加（n_sim = n_obs），同一归一化后投影进观测 PCA
  4. 双细胞得分 = 每个观测细胞 k = max(5, round(0.5·sqrt(n))) 近邻中模拟双细胞的比例
     （联合空间按观测轴标准化）
  5. 调用阈值 = 10x 期望双细胞率（≈0.8% × n/1000）对应得分分位 —— scDblFinder 同款做法；
     Otsu 仅作诊断字段留存（该数据得分分布近单峰，Otsu 会过度调用，不作推断口径）
  6. 分期比较：得分分布的 KS 检验（adenoma vs normal / adenoma vs carcinoma）
     —— 若腺瘤期得分不右移，即无"腺瘤双细胞更多"的证据
  7. 敏感性分解：剔除得分 top exp_rate（≈1%）与 top 5% 两档后，
     重算 sc161277_scan.py 三层读数（diff_abs 分化签名 T_all / U_epi / 上皮占比 P，
     另含 paneth 与 immune_U_epi 交叉污染校验），逐患者对比过滤前后

用法：python _supp_doublet.py   （纯本地 mtx.gz；约 1 分钟）
输出：_supp_doublet.json
"""
import gzip, io, json, os

import numpy as np
import scipy.io as sio
from scipy import stats as sps
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
D = os.path.join(HERE, "sc161277")

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
UMI_MIN, MT_MAX = 500, 0.20
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDH1"]
EPI_MIN_FRAC = 0.005
SIGS = {
    "diff_abs": ["ALPI", "FABP1", "FABP2", "CA1", "CA2", "GUCA2A", "GUCA2B", "CDX2", "CEACAM7"],
    "paneth":   ["DEFA5", "DEFA6", "LYZ", "REG3G", "PLA2G2A", "ITLN1"],
}
IMMUNE = ["PTPRC", "CD3D", "CD8A", "MS4A1", "LYZ", "LST1", "NKG7"]
N_HVG, N_PC, RNG_SEED = 3000, 50, 20261123

rng = np.random.default_rng(RNG_SEED)


def stem_of(gsm):
    for fn in os.listdir(D):
        if fn.startswith(gsm) and fn.endswith("_matrix.mtx.gz"):
            return fn[: -len("_matrix.mtx.gz")]
    return None


def otsu_threshold(scores):
    hist, edges = np.histogram(scores, bins=256)
    p = hist / hist.sum()
    omega = np.cumsum(p)
    mu = np.cumsum(p * np.arange(256))
    mu_t = mu[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b2 = (mu_t * omega - mu) ** 2 / (omega * (1 - omega))
    sigma_b2[~np.isfinite(sigma_b2)] = -1
    k = int(np.argmax(sigma_b2))
    return float((edges[k] + edges[k + 1]) / 2.0)


def process_sample(gsm, stem):
    feat = []
    with gzip.open(os.path.join(D, stem + "_features.tsv.gz"), "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            feat.append(p[1] if len(p) >= 2 else p[0])
    with gzip.open(os.path.join(D, stem + "_barcodes.tsv.gz"), "rt") as f:
        barcodes = [l.strip() for l in f if l.strip()]
    with gzip.open(os.path.join(D, stem + "_matrix.mtx.gz"), "rb") as f:
        X = sio.mmread(f)                       # genes × cells (COO)
    X = X.tocsc().astype(np.float32)            # 列 = 细胞

    sym = np.array(feat)
    ncell = X.shape[1]
    cell_total = np.asarray(X.sum(axis=0)).ravel()
    mt_mask = np.char.startswith(sym, "MT-")
    cell_mt = np.asarray(X[mt_mask].sum(axis=0)).ravel()
    epi_rows = np.where(np.isin(sym, EPI))[0]
    epi_umis = np.asarray(X[epi_rows].sum(axis=0)).ravel()

    keep = np.where((cell_total >= UMI_MIN) & (cell_mt / np.maximum(cell_total, 1) <= MT_MAX))[0]
    epi = keep[epi_umis[keep] / cell_total[keep] >= EPI_MIN_FRAC]
    non_epi = np.setdiff1d(keep, epi)

    # --- 基因过滤（在保留细胞中 ≥3 细胞表达） ---
    Xk = X[:, keep].tocsr()
    n_expr = (Xk > 0).sum(axis=1).A1
    gene_idx = np.where(n_expr >= 3)[0]
    Xg = Xk[gene_idx].tocsc()                   # genes_f × n_keep
    sym_f = sym[gene_idx]

    totals = np.asarray(Xg.sum(axis=0)).ravel()
    med_total = float(np.median(totals))
    scale = med_total / np.maximum(totals, 1.0)

    # --- HVG（log1p 后列方差 top N_HVG） ---
    Xl = Xg.tocsr().astype(np.float32)
    Xl.data = np.log1p(Xl.data * scale[Xl.nonzero()[1]])
    g_mean = np.asarray(Xl.mean(axis=1)).ravel()
    g_var = np.asarray(Xl.multiply(Xl).mean(axis=1)).ravel() - g_mean ** 2
    hvg = np.argsort(-g_var)[:N_HVG]

    Xh = np.asarray(Xl[hvg].T.todense(), dtype=np.float32)  # n_keep × N_HVG
    pca = PCA(n_components=min(N_PC, len(keep) - 1), svd_solver="full")
    Z_obs = pca.fit_transform(Xh)

    # --- 模拟双细胞：raw 相加 → 同一归一化 → HVG → PCA 投影 ---
    n_obs = len(keep)
    n_keep_n = Xg.shape[1]
    ii = rng.integers(0, n_keep_n, size=n_obs)
    jj = (ii + rng.integers(1, n_keep_n, size=n_obs)) % n_keep_n
    sim_raw = (Xg[:, ii] + Xg[:, jj]).tocsr().astype(np.float32)
    st = np.asarray(sim_raw.sum(axis=0)).ravel()
    sc_scale = med_total / np.maximum(st, 1.0)
    sim_raw.data = np.log1p(sim_raw.data * sc_scale[sim_raw.nonzero()[1]])
    Z_sim = pca.transform(np.asarray(sim_raw[hvg].T.todense(), dtype=np.float32))

    # --- kNN sim 占比得分（观测轴标准化后的联合空间） ---
    sd = Z_obs.std(axis=0)
    sd[sd == 0] = 1.0
    Zo = Z_obs / sd
    Zs = Z_sim / sd
    k = max(5, int(round(0.5 * np.sqrt(n_obs))))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(np.vstack([Zo, Zs]))
    _, idx = nn.kneighbors(Zo)                  # 第一列是自己
    score = (idx[:, 1:] >= n_obs).mean(axis=1)

    exp_rate = min(0.10, 0.008 * n_obs / 1000.0)     # 10x 期望双细胞率
    thr_quant = float(np.quantile(score, 1 - exp_rate))
    called = score >= thr_quant                       # 主口径：按期望率调用
    thr_otsu = otsu_threshold(score)                  # 诊断字段（不作推断口径）

    # --- 敏感性分解：剔除 called（≈期望率）与 top5% 两档，重算三层读数 ---
    sig_rows = {name: np.where(np.isin(sym, genes))[0] for name, genes in SIGS.items()}
    immune_rows = np.where(np.isin(sym, IMMUNE))[0]

    def sig_mean(cells, rows):
        if len(cells) == 0 or len(rows) == 0:
            return None
        sub = X[:, cells][rows, :].tocsc()
        t = cell_total[cells]                    # ★ 全基因 UMI 总量（与 sc161277_scan.py 一致）
        val = np.log1p(sub.multiply(1e4 / np.maximum(t, 1.0)))
        return float(np.asarray(val.mean(axis=0)).ravel().mean())

    def cell_mask_of(keep_arr, excluded_pos):
        """excluded_pos：keep 数组中被判为双细胞的行位置"""
        return keep_arr[~np.isin(np.arange(len(keep_arr)), excluded_pos)]

    excl_called = np.where(called)[0]
    excl_top5 = np.where(score >= np.quantile(score, 0.95))[0]

    out = {
        "gsm": gsm, "stage": SAMPLES[gsm][1], "patient": SAMPLES[gsm][0],
        "n_cells_raw": int(ncell), "n_cells_qc": int(len(keep)),
        "n_epi_qc": int(len(epi)),
        "knn_k": k, "n_sim": int(n_obs),
        "expected_rate": round(exp_rate, 4),
        "threshold_score": round(thr_quant, 4),
        "n_doublet_called": int(called.sum()),
        "doublet_rate_called": round(float(called.mean()), 4),
        "score_mean": round(float(score.mean()), 4),
        "score_median": round(float(np.median(score)), 4),
        "score_p95": round(float(np.percentile(score, 95)), 4),
        "otsu_threshold_diagnostic": round(thr_otsu, 4),
        "n_called_in_epi": int(called[np.isin(np.arange(len(keep)), [ {c: i for i, c in enumerate(keep)}[c] for c in epi ])].sum()) if len(epi) else 0,
    }

    decomp = {}
    for tag, excl in (("called", excl_called), ("top5", excl_top5)):
        keepnd = cell_mask_of(keep, excl)
        pos_of = {c: i for i, c in enumerate(keep)}
        epi_pos = np.array([pos_of[c] for c in epi])
        epi_nd = epi[~np.isin(epi_pos, excl)]
        non_nd = np.setdiff1d(keepnd, epi_nd)
        d = {}
        for name in SIGS:
            d[name] = {
                "T_all_raw": sig_mean(keep, sig_rows[name]),
                "T_all_filtered": sig_mean(keepnd, sig_rows[name]),
                "U_epi_raw": sig_mean(epi, sig_rows[name]),
                "U_epi_filtered": sig_mean(epi_nd, sig_rows[name]),
            }
        d["immune_U_epi_raw"] = sig_mean(epi, immune_rows)
        d["immune_U_epi_filtered"] = sig_mean(epi_nd, immune_rows)
        d["epi_frac_raw"] = float(len(epi)) / len(keep)
        d["epi_frac_filtered"] = float(len(epi_nd)) / len(keepnd) if len(keepnd) else None
        d["n_qc_after_filter"] = int(len(keepnd))
        decomp[tag] = d
    out["decomposition"] = decomp
    return out, score


results, all_scores = [], {}
for gsm, (pat, stage) in sorted(SAMPLES.items()):
    stem = stem_of(gsm)
    if stem is None:
        print("跳过（无矩阵）:", gsm, flush=True)
        continue
    print("处理", gsm, pat, stage, "…", flush=True)
    r, score = process_sample(gsm, stem)
    results.append(r)
    all_scores[gsm] = score.tolist()

# 与 scan JSON 的 QC 一致性断言
scan = {r["gsm"]: r for r in json.load(io.open("sc161277_scan.json", encoding="utf-8"))}
for r in results:
    ref = scan[r["gsm"]]
    assert r["n_cells_qc"] == ref["n_cells_qc"], ("QC 数不一致", r["gsm"])
    assert r["n_epi_qc"] == ref["_epi_idx_n"], ("上皮数不一致", r["gsm"])
print("QC 一致性断言：%d 样本全过" % len(results))

# 分期汇总：双细胞率（按期望率调用）+ 得分分布 KS 检验
summary = {}
for stage in ("normal", "adenoma", "carcinoma"):
    rs = [r for r in results if r["stage"] == stage]
    if rs:
        summary[stage] = {
            "n_samples": len(rs),
            "n_cells_qc_total": sum(r["n_cells_qc"] for r in rs),
            "n_doublet_called_total": sum(r["n_doublet_called"] for r in rs),
            "doublet_rate_called": round(sum(r["n_doublet_called"] for r in rs) /
                                         max(sum(r["n_cells_qc"] for r in rs), 1), 4),
            "score_median": round(float(np.median([x for r in rs
                                                   for x in all_scores[r["gsm"]]])), 4),
            "score_mean": round(float(np.mean([x for r in rs
                                               for x in all_scores[r["gsm"]]])), 4),
        }
sc_by_stage = {s: np.array([x for r in results if r["stage"] == s
                            for x in all_scores[r["gsm"]]])
               for s in ("normal", "adenoma", "carcinoma")}
ks_tests = {}
for a, b in (("adenoma", "normal"), ("adenoma", "carcinoma")):
    ks = sps.ks_2samp(sc_by_stage[a], sc_by_stage[b])
    ks_tests["%s_vs_%s" % (a, b)] = {
        "statistic": round(float(ks.statistic), 4),
        "p": float(ks.pvalue),
        "median_score_adenoma": summary["adenoma"]["score_median"],
        "median_score_other": summary[b]["score_median"],
    }
    print("KS %s vs %s: D=%.3f p=%.3g (median %.3f vs %.3f)"
          % (a, b, ks.statistic, ks.pvalue,
             summary["adenoma"]["score_median"], summary[b]["score_median"]))

print("\n分期双细胞率与得分：")
for s, v in summary.items():
    print("  %-10s n样本=%d QC细胞=%d 调用双细胞=%d (%.2f%%) score中位=%.3f"
          % (s, v["n_samples"], v["n_cells_qc_total"], v["n_doublet_called_total"],
             100 * v["doublet_rate_called"], v["score_median"]))

json.dump({"meta": {"method": "scrublet-style simulation + kNN sim-fraction score; threshold at "
                    "10x expected doublet rate (0.8% per 1000 cells); QC identical to "
                    "sc161277_scan.py; Otsu kept as diagnostic only (score distribution "
                    "near-unimodal, Otsu over-calls)",
                    "seed": RNG_SEED, "n_hvg": N_HVG, "n_pc": N_PC},
           "summary_by_stage": summary,
           "ks_tests": ks_tests,
           "samples": results,
           "scores_all_cells": all_scores},
          io.open("_supp_doublet.json", "w", encoding="utf-8"), ensure_ascii=False)
print("\n→ _supp_doublet.json")
