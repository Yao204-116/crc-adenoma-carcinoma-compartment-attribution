# -*- coding: utf-8 -*-
"""补充分析 ②：置换检验（B 类补分析，防"只报参数检验"指控）

三项内容：
  A. 逐基因方向一致率（上调分量 123/210 = 58.6%；下调分量 96/103 = 93.2%）
       口径 = manuscript_facts.py 5b 节：x = mean(d_GSE4183, d_GSE117606)，y = d_GSE15960_LCM，
       同向 = x*y > 0（严格，0 不计同向）。
       检验：① 精确二项（复核摘要 P=0.016）；② 符号翻转置换 20k；
             ③ 边际保持的基因配对打乱 20k（把 y 向量整体置换，保留 LCM 上调基因占比）——
               这是更保守的零假设（LCM 侧 58.6% 基因上调的边际不平衡会让 ≥50% 的朴素零假设失真）。
  B. 配对计数（全签名 53/56；上皮分化·核心8 面板 8/56）
       口径 = gse117606_paired.py 分析 B：Δ = 癌 − 腺瘤（同患者），单位=队列内 z。
       检验：逐对符号翻转置换 20k（等价于精确符号检验的 MC 版）+ 双侧 P。
       点位断言：必须复现 gse117606_paired_results.json 的 n_pair/n_pos。
用法：python _supp_permutation.py   （离线缓存；载入 GSE117606 约需 1-2 分钟）
输出：_supp_permutation.json
"""
import csv as _csv
import io, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
rng = np.random.default_rng(20261123)
B = 20000

# ===========================================================================
# A. 逐基因方向一致率（附表 1 驱动，无需重载表达矩阵）
# ===========================================================================
rows = list(_csv.DictReader(io.open("suppl_table1_313genes.csv", encoding="utf-8-sig")))
comp_data = {}
for comp, key in ((1, "UP"), (-1, "DN")):
    xs, ys = [], []
    for r in rows:
        if (1 if float(r["mean_logFC"]) > 0 else -1) != comp:
            continue
        # 5b 精确口径：全组织 delta 至少一个存在即入组（均值取可得的）
        t = [float(r[c]) for c in ("d_GSE4183", "d_GSE117606")
             if r.get(c) not in ("", None, "None")]
        y = r.get("d_GSE15960_LCM")
        if not t or y in ("", None, "None"):
            continue
        xs.append(sum(t) / len(t))
        ys.append(float(y))
    comp_data[key] = (np.array(xs), np.array(ys))

from math import comb

def binom_two_sided(k, n, p=0.5):
    def pmf(x):
        return comb(n, x) * p ** x * (1 - p) ** (n - x)
    obs = pmf(k)
    return sum(pmf(x) for x in range(n + 1) if pmf(x) <= obs + 1e-12)

# ===========================================================================
# A-③ 全基因组基线（GSE15960 内）：123/210 是否高于「任意基因在 LCM 中上调」的基础率
#   零假设：UP 分量基因在 LCM 中的上调率 = 全基因组上调率 p0（同口径 drift 校正后 Δ）
#   注：边际保持打乱检验循环（x 全为正 → 零假设均值=观测值），仅留档说明，不作推断。
# ===========================================================================
_src = io.open("manuscript_facts.py", encoding="utf-8").read()
_ns2 = {"__name__": "_mfdefs", "__file__": os.path.join(HERE, "manuscript_facts.py")}
exec(_src[:_src.index('print("载入队列")')], _ns2)
_c15960 = _ns2["load_cohort"]("GSE15960")
zmap = _c15960["zmap_adj"]
lab = _c15960["labels"]
ia = [i for i, l in enumerate(lab) if l == "adenoma"]
ic = [i for i, l in enumerate(lab) if l == "carcinoma"]
import statistics as _st
ds = []
for g, z in zmap.items():
    ca = [z[i] for i in ic if i < len(z)]
    ad = [z[i] for i in ia if i < len(z)]
    if len(ca) >= 2 and len(ad) >= 2:
        ds.append(_st.mean(ca) - _st.mean(ad))
p0 = sum(1 for v in ds if v > 0) / len(ds)
print("A-③ GSE15960 全基因组 Δ(癌−腺瘤)>0 的基础率 p0 = %.3f（n=%d 基因）" % (p0, len(ds)), flush=True)

perm_out = {"gene_concordance": {}, "paired_counts": {},
            "genome_wide_baseline": {"p0": round(p0, 4), "n_genes": len(ds)}}
EXPECT = {"UP": (210, 123), "DN": (103, 96)}
for key, (x, y) in comp_data.items():
    n = len(x)
    obs = int(np.sum(x * y > 0))
    exp_n, exp_k = EXPECT[key]
    assert (n, obs) == (exp_n, exp_k), ("一致率真值断言失败", key, n, obs, exp_n, exp_k)
    # ① 精确二项（vs 0.5）
    p_binom = binom_two_sided(obs, n)
    # ② 符号翻转置换
    flips = (rng.random((B, n)) < 0.5)
    cnt_flip = np.zeros(B, dtype=int)
    ytile = np.tile(y, (B, 1))
    yflip = np.where(flips, -ytile, ytile)
    cnt_flip = np.sum((x[None, :] * yflip) > 0, axis=1)
    mean_null = cnt_flip.mean()
    p_flip = float(np.mean(np.abs(cnt_flip - mean_null) >= abs(obs - mean_null)))
    # ③ 边际保持打乱（置换 y）
    cnt_shuf = np.empty(B, dtype=int)
    for b in range(B):
        yy = y[rng.permutation(n)]
        cnt_shuf[b] = int(np.sum(x * yy > 0))
    p_shuf = float(np.mean(np.abs(cnt_shuf - cnt_shuf.mean()) >= abs(obs - cnt_shuf.mean())))
    perm_out["gene_concordance"][key] = {
        "n_genes": n, "observed_concordant": obs,
        "pct": round(100.0 * obs / n, 1),
        "binomial_two_sided_p_vs_0.5": p_binom,
        "binomial_two_sided_p_vs_p0": binom_two_sided(obs, n, p0),
        "signflip_permutation": {"B": B, "null_mean": round(float(mean_null), 2),
                                 "two_sided_p": p_flip},
        "marginal_shuffle_permutation": {"B": B, "null_mean": round(float(cnt_shuf.mean()), 2),
                                         "two_sided_p": p_shuf,
                                         "note": "circular: tissue-side x is all-up (UP) / all-but-4-down "
                                                 "(DN), so the shuffle null reproduces the y-marginal and "
                                                 "the observed value; kept for transparency only"},
    }
    print("A. %s 分量：%d/%d = %.1f%%" % (key, obs, n, 100 * obs / n))
    print("   二项双侧 p = %.4f | 符号翻转 p = %.4f（零假设均值 %.1f）| 打乱 p = %.4f（零假设均值 %.1f）"
          % (p_binom, p_flip, mean_null, p_shuf, cnt_shuf.mean()), flush=True)

# ===========================================================================
# B. 配对计数（重用 gse117606_paired.py 的构造段）
# ===========================================================================
src = io.open("gse117606_paired.py", encoding="utf-8").read()
ns = {"__name__": "_gpdefs", "__file__": os.path.join(HERE, "gse117606_paired.py")}
exec(src[:src.index("def auc_wrap(")], ns)
score_sig, score = ns["score_sig"], ns["score"]
pairs_ac, byp = ns["pairs_ac"], ns["byp"]
ALLPANELS = ns["ALLPANELS"]

saved = json.load(io.open("gse117606_paired_results.json", encoding="utf-8"))

def pair_deltas(fn):
    out = []
    for p in pairs_ac:
        sa, sc = fn(byp[p]["adenoma"][0]), fn(byp[p]["carcinoma"][0])
        if sa is None or sc is None:
            continue
        out.append(sc - sa)
    return np.array(out)

for tag, fn, json_key, expect_pos in [
    ("full_signature", score_sig, ("B_paired", "全签名"), 53),
    ("epi_diff_core8", lambda i: score(i, ALLPANELS["上皮分化·核心8"]),
     ("B_panel_paired", "上皮分化·核心8"), 8),
    ("CAF_panel", lambda i: score(i, ALLPANELS["间质·激活/CAF"]),
     ("B_panel_paired", "间质·激活/CAF"), 53),
]:
    d = pair_deltas(fn)
    n = len(d)
    npos = int(np.sum(d > 0))
    ref = saved[json_key[0]][json_key[1]]
    n_ref = ref.get("n_pair", ref.get("n"))
    assert (n_ref, ref["n_pos"]) == (n, npos), ("配对真值断言失败", tag, n, npos, ref)
    # 符号翻转置换（双侧，以零假设均值 n/2 为中心）
    sims = np.abs(rng.binomial(n, 0.5, B) - n / 2.0)
    p_perm = float(np.mean(sims >= abs(npos - n / 2.0)))
    perm_out["paired_counts"][tag] = {
        "n_pairs": n, "n_pos": npos,
        "sign_test_p_saved": ref.get("sign_test_p"),
        "wilcoxon_p_saved": ref.get("wilcoxon_p"),
        "signflip_permutation_two_sided_p": p_perm,
        "median_delta": float(np.median(d)),
    }
    print("B. %s：%d/%d 癌>腺瘤 | 置换双侧 p = %.5f（稿件符号检验 p = %s）"
          % (tag, npos, n, p_perm, ref.get("sign_test_p")), flush=True)

json.dump(perm_out, io.open("_supp_permutation.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n→ _supp_permutation.json")
