# -*- coding: utf-8 -*-
"""补充分析 ①：关键 AUC 的 bootstrap 95% CI（B 类补分析，防审稿人"无区间"指控）

口径 = 稿件 Methods 无向打分：
  - 逐基因 z（队内全体样本标准化）→ 逐样本扣全局水平（zmap_adj，漂移校正）
  - 基因集得分 = 成员基因 z 的均值（无向，与 final_panel_matrix.csv 一致）
  - AUC = P(score | carcinoma > score | adenoma)，方向固定，不做 max(v,1-v) 翻转
  - 分层自助：腺瘤、癌组内分别有放回重抽样，B=3000，百分位法 95% CI
复现性：必须逐位复现 final_panel_matrix.csv 的点位（脚本内置断言）。
用法：python _supp_bootstrap_auc_ci.py   （离线，只用 cohort_hunt/geo_cache）
输出：_supp_bootstrap_auc_ci.json
"""
import io, json, os, re, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

# ---- 复用 manuscript_facts.py 的加载段（截到"载入队列"横幅之前，纯定义无副作用）----
_src = io.open("manuscript_facts.py", encoding="utf-8").read()
_cut = _src.index('print("载入队列")')
_ns = {"__name__": "_mfdefs", "__file__": os.path.join(HERE, "manuscript_facts.py")}
exec(_src[:_cut], _ns)
load_cohort = _ns["load_cohort"]
COHORTS = _ns["COHORTS"]
COHORT_NOTE = _ns["COHORT_NOTE"]
PANELS = _ns["PANELS"]
auc_ref = _ns["auc"]

import csv as _csv
SIGINFO = {r["gene"].upper(): r for r in
           _csv.DictReader(io.open("consensus_signature_tier.csv", encoding="utf-8-sig"))}
UP = [g for g, r in SIGINFO.items() if float(r["mean_logFC"]) > 0]
DN = [g for g, r in SIGINFO.items() if float(r["mean_logFC"]) < 0]
SETS = {
    "signature_313": ("Full 313-gene signature", sorted(SIGINFO)),
    "UP_component": ("Up-regulated component (210)", sorted(UP)),
    "DN_component": ("Down-regulated component (103)", sorted(DN)),
    "pure_epi_baseline7": ("Pure-epithelial baseline (7)", sorted(PANELS["纯上皮基线·7"])),
    "epi_diff_core8": ("Epithelial differentiation core (8)", sorted(PANELS["上皮分化·核心8"])),
}

B = 3000
rng = np.random.default_rng(20261123)


def auc_fast(c_scores, a_scores):
    """Mann-Whitney AUC = P(c > a) + 0.5 P(c == a)，平均秩实现（与 sig_validate.auc 等价）"""
    c = np.asarray(c_scores, dtype=float)
    a = np.asarray(a_scores, dtype=float)
    n_c, n_a = len(c), len(a)
    comb = np.concatenate([c, a])
    order = comb.argsort(kind="mergesort")
    ranks = np.empty(len(comb))
    ranks[order] = np.arange(1, len(comb) + 1)
    # 并列取平均秩
    sc = comb[order]
    i = 0
    while i < len(sc):
        j = i
        while j + 1 < len(sc) and sc[j + 1] == sc[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return (ranks[:n_c].sum() - n_c * (n_c + 1) / 2.0) / (n_c * n_a)


OUT = {"meta": {"B": B, "seed": 20261123, "caliber": "undirected mean-z, drift-corrected (zmap_adj); "
                "AUC = P(carcinoma > adenoma), orientation fixed by the data (no max-flip)",
                "note": "z / drift-correction computed once on the full cohort (incl. normals); "
                        "bootstrap resamples adenoma and carcinoma specimens with replacement within group"},
       "results": []}

print("载入队列（离线缓存）…", flush=True)
C = {}
for gse in COHORTS:
    c = load_cohort(gse)
    assert c is not None, gse + " 载入失败"
    C[gse] = c
    from collections import Counter
    cnt = Counter(c["labels"])
    print("  %-11s 基因 %5d | 腺瘤 %d 癌 %d 正常 %d" % (gse, len(c["zmap"]),
          cnt.get("adenoma", 0), cnt.get("carcinoma", 0), cnt.get("normal", 0)), flush=True)

for gse in COHORTS:
    zmap = C[gse]["zmap_adj"]
    lab = C[gse]["labels"]
    ia = np.array([i for i, l in enumerate(lab) if l == "adenoma"])
    ic = np.array([i for i, l in enumerate(lab) if l == "carcinoma"])
    for key, (_desc, genes) in SETS.items():
        gs = [g for g in genes if g in zmap]
        if len(gs) < 3:
            continue
        # 逐样本得分（一次算好；bootstrap 只重排）
        M = np.array([[zmap[g][i] if i < len(zmap[g]) else 0.0 for g in gs]
                      for i in range(len(lab))], dtype=float)
        s_all = M.mean(axis=1)
        s_a, s_c = s_all[ia], s_all[ic]
        pt = auc_fast(s_c, s_a)
        # 复现断言（容差 1e-9）
        ref = auc_ref(list(s_c), list(s_a))
        assert abs(pt - ref) < 1e-9, (gse, key, pt, ref)
        # 分层 bootstrap
        n_a, n_c_ = len(ia), len(ic)
        boots = np.empty(B)
        for b in range(B):
            ba = rng.choice(n_a, n_a, replace=True)
            bc = rng.choice(n_c_, n_c_, replace=True)
            boots[b] = auc_fast(s_c[bc], s_a[ba])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        OUT["results"].append({
            "cohort": gse, "kind": COHORT_NOTE[gse][0], "gene_set": key,
            "n_adenoma": n_a, "n_carcinoma": n_c_, "n_genes": len(gs),
            "auc_point": round(float(pt), 4),
            "ci95_low": round(float(lo), 4), "ci95_high": round(float(hi), 4),
            "ci_contains_0.5": bool(lo <= 0.5 <= hi),
        })
        print("  %-11s %-22s AUC %.3f  [%.3f, %.3f]%s"
              % (gse, key, pt, lo, hi, "  ←含0.5" if lo <= 0.5 <= hi else ""), flush=True)

json.dump(OUT, io.open("_supp_bootstrap_auc_ci.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n→ _supp_bootstrap_auc_ci.json（%d 条）" % len(OUT["results"]))
