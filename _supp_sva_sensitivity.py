# -*- coding: utf-8 -*-
"""补充分析 ④：SVA 式敏感性分析 + 敏感性图（B 类补分析）

目的：排除"替代变量（批次/混杂）驱动签名"的指控。

方法（SVA-lite，逐队列独立）：
  1. 从 zmap_adj（逐基因 z + 全局漂移校正，与稿件口径一致）构建 基因×样本 矩阵
  2. 逐基因回归掉阶段效应（normal/adenoma/carcinoma/ssa/other 单热），取残差
  3. 残差 SVD → 前 5 个替代变量 SV（按构造与阶段正交，捕捉批次结构）
  4. 敏感性：k = 0..5，把前 k 个 SV（+截距）从原始矩阵逐基因回归掉，
     重算：① 313 签名无向 AUC（癌 vs 腺瘤）② 上调分量保留率 ③ 下调分量放大幅度
           ④ 逐基因方向一致率（58.6% / 93.2% 口径 = manuscript_facts 5b）
  5. 每个 SV 与阶段的关联强度（one-way ANOVA η²）——佐证 SV 抓的是批次不是阶段

用法：python _supp_sva_sensitivity.py   （离线缓存；约 2-3 分钟）
输出：_supp_sva_sensitivity.json + figS_sva_sensitivity.png（300 dpi，英文投稿图）
"""
import csv as _csv
import io, json, os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

# ---- 复用 manuscript_facts.py 的加载段（纯定义） ----
_src = io.open("manuscript_facts.py", encoding="utf-8").read()
_ns = {"__name__": "_mfdefs", "__file__": os.path.join(HERE, "manuscript_facts.py")}
exec(_src[:_src.index('print("载入队列")')], _ns)
load_cohort = _ns["load_cohort"]
COHORTS = _ns["COHORTS"]
PANELS = _ns["PANELS"]

SIGINFO = {r["gene"].upper(): r for r in
           _csv.DictReader(io.open("consensus_signature_tier.csv", encoding="utf-8-sig"))}
UP = sorted(g for g, r in SIGINFO.items() if float(r["mean_logFC"]) > 0)
DN = sorted(g for g, r in SIGINFO.items() if float(r["mean_logFC"]) < 0)
SIG = UP + DN
EPI7 = sorted(PANELS["纯上皮基线·7"])
K_MAX = 5

print("载入队列（离线缓存）…", flush=True)
C = {}
for gse in COHORTS:
    C[gse] = load_cohort(gse)
    print("  %s ok" % gse, flush=True)


def build_matrix(gse):
    c = C[gse]
    genes = sorted(c["zmap_adj"].keys())
    nsamp = len(c["labels"])
    Z = np.array([c["zmap_adj"][g] for g in genes], dtype=np.float64)
    return genes, Z, c["labels"]


def sv_extract(Z, labels, k_max=5):
    """阶段残差的 SVD → 样本荷载 SV（正交于阶段）"""
    cats = sorted(set(labels))
    D = np.zeros((len(labels), len(cats)))
    for i, l in enumerate(labels):
        D[i, cats.index(l)] = 1.0
    # 逐基因回归掉阶段 → 残差
    beta, *_ = np.linalg.lstsq(D, Z.T, rcond=None)
    Res = (Z.T - D @ beta).T                              # genes × samples
    M = Res.T                                             # samples × genes
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    SVs = U[:, :k_max] * S[:k_max]                        # 样本荷载
    # SV 与阶段的关联（one-way ANOVA eta²，仅三类主阶段）
    keep = [i for i, l in enumerate(labels) if l in ("normal", "adenoma", "carcinoma")]
    assoc = []
    for j in range(k_max):
        x = SVs[:, j]
        groups = [x[[i for i in keep if labels[i] == s]] for s in ("normal", "adenoma", "carcinoma")]
        grand = x[keep].mean()
        ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups if len(g))
        ss_t = ((x[keep] - grand) ** 2).sum()
        assoc.append(float(ss_b / ss_t) if ss_t > 0 else 0.0)
    return SVs, assoc


def residualize(Z, SVs, k):
    if k == 0:
        return Z
    D = np.hstack([np.ones((SVs.shape[0], 1)), SVs[:, :k]])
    beta, *_ = np.linalg.lstsq(D, Z.T, rcond=None)
    return (Z.T - D @ beta).T


def auc_fast(c_scores, a_scores):
    comb = np.concatenate([c_scores, a_scores])
    order = comb.argsort(kind="mergesort")
    ranks = np.empty(len(comb))
    ranks[order] = np.arange(1, len(comb) + 1)
    sc = comb[order]
    i = 0
    while i < len(sc):
        j = i
        while j + 1 < len(sc) and sc[j + 1] == sc[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    n_c = len(c_scores)
    return (ranks[:n_c].sum() - n_c * (n_c + 1) / 2.0) / (n_c * len(a_scores))


KS = {k: 0 for k in range(K_MAX + 1)}   # k=0 复现校验
OUT = {"meta": {"method": "SVA-lite: per-cohort stage-residual SVD -> top-5 surrogate variables; "
                "recompute after regressing out first k SVs (+intercept); downstream caliber "
                "identical to manuscript (undirected mean-z, drift-corrected)"},
       "sv_stage_assoc": {}, "sensitivity": {}}

for gse in COHORTS:
    genes, Z, labels = build_matrix(gse)
    gidx = {g: i for i, g in enumerate(genes)}
    SVs, assoc = sv_extract(Z, labels, K_MAX)
    OUT["sv_stage_assoc"][gse] = {"eta2_sv1_sv5": [round(a, 4) for a in assoc]}
    ia = np.array([i for i, l in enumerate(labels) if l == "adenoma"])
    ic = np.array([i for i, l in enumerate(labels) if l == "carcinoma"])

    sig_rows = [gidx[g] for g in SIG if g in gidx]
    up_rows = [gidx[g] for g in UP if g in gidx]
    dn_rows = [gidx[g] for g in DN if g in gidx]

    per_k = {}
    for k in range(K_MAX + 1):
        Zk = residualize(Z, SVs, k)
        # 签名 AUC
        s = Zk[sig_rows].mean(axis=0)
        auc = auc_fast(s[ic], s[ia])
        # 逐基因 delta（5b 口径）
        def delta(rows):
            return Zk[rows][:, ic].mean(axis=1) - Zk[rows][:, ia].mean(axis=1)
        per_k[k] = {"auc_313": round(float(auc), 4),
                    "d_up_mean": float(delta(up_rows).mean()),
                    "d_dn_mean": float(delta(dn_rows).mean()),
                    "up_frac_pos": round(float((delta(up_rows) > 0).mean()), 4),
                    "dn_frac_pos": round(float((delta(dn_rows) > 0).mean()), 4)}
    OUT["sensitivity"][gse] = per_k
    print("%s: AUC(313) k=0..5 = %s" % (gse, [per_k[k]["auc_313"] for k in range(K_MAX + 1)]),
          flush=True)

# ---- 复现校验：k=0 时 AUC 应与 final_panel_matrix.csv 一致 ----
ref = {r["基因集"]: r for r in _csv.DictReader(io.open("final_panel_matrix.csv", encoding="utf-8-sig"))}
for gse in COHORTS:
    got = OUT["sensitivity"][gse][0]["auc_313"]
    exp = float(ref["313签名"]["%s_auc" % gse])
    assert abs(got - exp) < 5e-5, ("k=0 复现失败", gse, got, exp)  # 存储已舍入到 4 位
print("k=0 复现断言（3 队列 AUC = final_panel_matrix.csv）：全过")

# ---- 跨队列汇总量（敏感性图的主对象） ----
def retention(k):
    """UP 保留率 % 与 DN 放大倍数 + 一致率（5b 口径，x=两全组织均值, y=LCM）"""
    rows = list(_csv.DictReader(io.open("suppl_table1_313genes.csv", encoding="utf-8-sig")))
    # 注意：附表的 d_ 值是 k=0 的；敏感性的逐基因 delta 需要重算 → 直接从 Z 重算
    return None

# 5b 口径敏感性：x = mean(d_4183, d_117606 可得者), y = d_GSE15960；同向 = x*y>0
comp_data = {}
for g, r in SIGINFO.items():
    comp_data[r["gene"].upper()] = r

Zs = {}
for gse in COHORTS:
    genes, Z, labels = build_matrix(gse)
    Zs[gse] = (genes, {g: i for i, g in enumerate(genes)}, Z, labels)

sv5 = {gse: sv_extract(Zs[gse][2], Zs[gse][3], K_MAX)[0] for gse in COHORTS}

# ---- 激进版上界：SV 直接从原始矩阵取（含阶段信号），仅存档，不作主口径 ----
def sv_extract_raw(Z, k_max=5):
    """不做阶段残差化：SVD 直接取原始矩阵的 top SV（可能携带阶段信号）"""
    M = Z.T
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    return U[:, :k_max] * S[:k_max]

def stage_eta2(sv, labels):
    keep = [i for i, l in enumerate(labels) if l in ("normal", "adenoma", "carcinoma")]
    x = sv
    grand = x[keep].mean()
    groups = [x[[i for i in keep if labels[i] == s]] for s in ("normal", "adenoma", "carcinoma")]
    ss_b = sum(len(g) * (g.mean() - grand) ** 2 for g in groups if len(g))
    ss_t = ((x[keep] - grand) ** 2).sum()
    return float(ss_b / ss_t) if ss_t > 0 else 0.0

sv_raw = {gse: sv_extract_raw(Zs[gse][2], K_MAX) for gse in COHORTS}
OUT["aggressive_overcorrection"] = {"sv_stage_eta2": {}, "auc_313": {}, "note":
    "SVs taken from the UNADJUSTED matrix (may carry stage signal). Removing them is an "
    "aggressive over-correction upper bound, not a valid null — the leading unsupervised "
    "directions largely ARE the biological contrast. Stored for transparency only."}
for gse in COHORTS:
    genes, gidx, Z, labels = Zs[gse]
    OUT["aggressive_overcorrection"]["sv_stage_eta2"][gse] = [
        round(stage_eta2(sv_raw[gse][:, j], labels), 3) for j in range(K_MAX)]
    ia = np.array([i for i, l in enumerate(labels) if l == "adenoma"])
    ic = np.array([i for i, l in enumerate(labels) if l == "carcinoma"])
    sig_rows = [gidx[g] for g in SIG if g in gidx]
    tmp = {}
    for k in range(K_MAX + 1):
        Zk = residualize(Z, sv_raw[gse], k)
        s = Zk[sig_rows].mean(axis=0)
        tmp[k] = round(float(auc_fast(s[ic], s[ia])), 4)
    OUT["aggressive_overcorrection"]["auc_313"][gse] = tmp
    print("激进版 %s: AUC k=0..5 = %s  (raw-SV 阶段 η² = %s)"
          % (gse, [tmp[k] for k in range(K_MAX + 1)],
             OUT["aggressive_overcorrection"]["sv_stage_eta2"][gse]), flush=True)

sens_5b = {}
for k in range(K_MAX + 1):
    # 逐队列残差化后的 delta
    dmap = {}
    for gse in COHORTS:
        genes, gidx, Z, labels = Zs[gse]
        Zk = residualize(Z, sv5[gse], k)
        ia = [i for i, l in enumerate(labels) if l == "adenoma"]
        ic = [i for i, l in enumerate(labels) if l == "carcinoma"]
        rows_up = [gidx[g] for g in UP if g in gidx]
        rows_dn = [gidx[g] for g in DN if g in gidx]
        dmap[gse] = {"UP": Zk[rows_up][:, ic].mean(axis=1) - Zk[rows_up][:, ia].mean(axis=1),
                     "DN": Zk[rows_dn][:, ic].mean(axis=1) - Zk[rows_dn][:, ia].mean(axis=1),
                     "up_genes": [g for g in UP if g in gidx],
                     "dn_genes": [g for g in DN if g in gidx]}
    # 5b 跨队列
    up_x, up_y, dn_x, dn_y = [], [], [], []
    for r in _csv.DictReader(io.open("suppl_table1_313genes.csv", encoding="utf-8-sig")):
        g = r["gene"]
        comp = 1 if float(r["mean_logFC"]) > 0 else -1
        key = "UP" if comp == 1 else "DN"
        if g not in dmap["GSE4183"]["up_genes"] and g not in dmap["GSE4183"]["dn_genes"]:
            continue
        # 在 k 残差下重算的 delta
        def dv(gse):
            lst, di = ("up_genes", "UP") if key == "UP" else ("dn_genes", "DN")
            if g not in dmap[gse][lst]:
                return None
            idx = dmap[gse][lst].index(g)
            return float(dmap[gse][key][idx])
        t = [dv("GSE4183"), dv("GSE117606")]
        t = [v for v in t if v is not None]
        y = dv("GSE15960")
        if not t or y is None:
            continue
        x = sum(t) / len(t)
        (up_x if comp == 1 else dn_x).append(x)
        (up_y if comp == 1 else dn_y).append(y)
    up_ret = float(np.mean(up_y) / np.mean(up_x)) if np.mean(up_x) else None
    dn_amp = float(np.mean(dn_y) / np.mean(dn_x))
    up_conc = float(np.mean([1.0 if x * y > 0 else 0.0 for x, y in zip(up_x, up_y)]))
    dn_conc = float(np.mean([1.0 if x * y > 0 else 0.0 for x, y in zip(dn_x, dn_y)]))
    sens_5b[k] = {"n_up": len(up_x), "n_dn": len(dn_x),
                  "up_retention": round(up_ret, 4), "dn_amplification": round(dn_amp, 4),
                  "up_concordance": round(up_conc, 4), "dn_concordance": round(dn_conc, 4)}
    print("k=%d: UP保留 %.3f | DN放大 %.3f | 一致率 UP %.3f / DN %.3f"
          % (k, up_ret, dn_amp, up_conc, dn_conc), flush=True)

# k=0 校验：31.6% / 2.04× / 58.6% / 93.2%
assert abs(sens_5b[0]["up_retention"] - 0.3158) < 0.01, sens_5b[0]
assert abs(sens_5b[0]["dn_amplification"] - 2.0376) < 0.02, sens_5b[0]
assert abs(sens_5b[0]["up_concordance"] - 0.5857) < 0.01, sens_5b[0]
assert abs(sens_5b[0]["dn_concordance"] - 0.9320) < 0.01, sens_5b[0]
print("k=0 复现断言（31.6% / 2.04× / 58.6% / 93.2%）：全过")

OUT["sensitivity_5b"] = sens_5b
json.dump(OUT, io.open("_supp_sva_sensitivity.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("→ _supp_sva_sensitivity.json")

# ================= 敏感性图（英文投稿图，300 dpi） =================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ks = list(range(K_MAX + 1))
fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))
colors = {"GSE4183": "#1a56b0", "GSE117606": "#c0392b", "GSE15960": "#2e7d32"}

ax = axes[0][0]
for gse in COHORTS:
    ys = [OUT["sensitivity"][gse][k]["auc_313"] for k in ks]
    ax.plot(ks, ys, "o-", color=colors[gse], label=gse)
ax.axhline(0.5, color="#888", ls="--", lw=1)
ax.set_xlabel("Number of surrogate variables removed (k)")
ax.set_ylabel("Signature AUC (carcinoma vs adenoma)")
ax.set_title("(a) Full-signature AUC is stable to surrogate-variable removal", fontsize=10)
ax.legend(fontsize=8, frameon=False)
ax.set_ylim(-0.05, 1.05)

ax = axes[0][1]
ax.plot(ks, [100 * sens_5b[k]["up_retention"] for k in ks], "o-", color="#1a56b0",
        label="Up component retained (%)")
ax.axhline(100, color="#888", ls="--", lw=1)
ax.set_xlabel("Number of surrogate variables removed (k)")
ax.set_ylabel("Retention of whole-tissue amplitude (%)")
ax.set_title("(b) Up-component amplitude retention (~32%)", fontsize=10)
ax.legend(fontsize=8, frameon=False)

ax = axes[1][0]
ax.plot(ks, [sens_5b[k]["dn_amplification"] for k in ks], "o-", color="#c0392b",
        label="Down component amplification (×)")
ax.axhline(1.0, color="#888", ls="--", lw=1)
ax.set_xlabel("Number of surrogate variables removed (k)")
ax.set_ylabel("Amplification relative to whole tissue (×)")
ax.set_title("(c) Down-component amplification (~2×)", fontsize=10)
ax.legend(fontsize=8, frameon=False)

ax = axes[1][1]
ax.plot(ks, [100 * sens_5b[k]["up_concordance"] for k in ks], "o-", color="#1a56b0",
        label="Up component")
ax.plot(ks, [100 * sens_5b[k]["dn_concordance"] for k in ks], "o-", color="#7b3fa0",
        label="Down component")
ax.axhline(50, color="#888", ls="--", lw=1)
ax.set_xlabel("Number of surrogate variables removed (k)")
ax.set_ylabel("Gene-level direction concordance (%)")
ax.set_title("(d) Direction concordance between whole tissue and LCM", fontsize=10)
ax.legend(fontsize=8, frameon=False)

for a in axes.ravel():
    a.set_xticks(ks)
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)

fig.suptitle("Sensitivity to hidden structure: top surrogate variables removed one at a time",
             fontsize=12, y=0.99)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig("figS_sva_sensitivity.png", dpi=300)
print("→ figS_sva_sensitivity.png")
