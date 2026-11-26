# -*- coding: utf-8 -*-
"""
R2 审稿 §3.2（M6 残留）：基因共表达导致方向一致率检验的 p 值高估？
按审稿人字面建议做两件事：
  ① 模块分块符号翻转置换（保留 WGCNA 共表达模块结构，块整体同翻）→ 重算 58.6% / 93.2% 的 p
  ② CAMERA 式方差膨胀：估基因间平均相关 ρ̄ → VIF=1+(n-1)ρ̄ → 二项 z 除以 sqrt(VIF)
口径与 _supp_permutation.py 完全一致（附表 1 驱动，5b "至少一个全组织 delta 存在"）。
断言：123/210、96/103 逐位复现。
"""
import io, os, csv, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

rows = list(csv.DictReader(io.open("suppl_table1_313genes.csv", encoding="utf-8-sig")))

def load_comp(comp):
    """comp=+1 → UP, −1 → DN；返回 genes, x(全组织均值 delta), y(LCM delta)"""
    genes, xs, ys = [], [], []
    for r in rows:
        if (1 if float(r["mean_logFC"]) > 0 else -1) != comp:
            continue
        t = [float(r[c]) for c in ("d_GSE4183", "d_GSE117606")
             if r.get(c) not in ("", None, "None")]
        y = r.get("d_GSE15960_LCM")
        if not t or y in ("", None, "None"):
            continue
        genes.append(r["gene"].upper())
        xs.append(sum(t) / len(t))
        ys.append(float(y))
    return genes, np.array(xs), np.array(ys)

UP_g, UP_x, UP_y = load_comp(1)
DN_g, DN_x, DN_y = load_comp(-1)
n_up, k_up = len(UP_g), int((UP_x * UP_y > 0).sum())
n_dn, k_dn = len(DN_g), int((DN_x * DN_y > 0).sum())
assert (n_up, k_up) == (210, 123), (n_up, k_up)
assert (n_dn, k_dn) == (103, 96), (n_dn, k_dn)
print("锚定 OK：UP %d/%d (58.6%%) | DN %d/%d (93.2%%)" % (k_up, n_up, k_dn, n_dn), flush=True)

# ---------------------------------------------------------------------------
# ① 模块分块符号翻转
# ---------------------------------------------------------------------------
mod = {}
_UPR = os.environ.get("AACA_UPSTREAM_R") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "upstream_r")
_WGC_DIR = os.path.join(_UPR, "02.WGCNA")
if not os.path.isfile(os.path.join(_WGC_DIR, "gene_module_assignment.csv")):
    _WGC_DIR = r"D:/2/02.WGCNA"  # developer-machine fallback (matches README §3)
with io.open(os.path.join(_WGC_DIR, "gene_module_assignment.csv"), encoding="utf-8-sig") as fh:
    rd = csv.DictReader(fh)
    cols = rd.fieldnames
    print("module csv 字段:", cols, flush=True)
    for r in rd:
        g = str(r["Gene"]).upper()
        mod[g] = r["Module"]
print("模块表覆盖 %d 基因" % len(mod), flush=True)

def block_labels(genes):
    """每基因的块 id：WGCNA 模块；无模块信息 → 单基因块"""
    lab, nxt = {}, 0
    for g in genes:
        m = mod.get(g)
        if m:
            lab[g] = "mod:" + m
        else:
            lab[g] = "solo:%d" % nxt
            nxt += 1
    return lab

def block_flip_pvalue(genes, xs, ys, B=20000, seed=20261125):
    lab = block_labels(genes)
    blocks = {}
    for g in genes:
        blocks.setdefault(lab[g], []).append(g)
    sizes = sorted((len(v) for v in blocks.values()), reverse=True)
    obs = int((xs * ys > 0).sum())
    n = len(genes)
    idx = {g: i for i, g in enumerate(genes)}
    block_idx = [[idx[g] for g in v] for v in blocks.values()]
    rng = np.random.default_rng(seed)
    ge = 0
    flip = np.ones(n)
    for b in range(B):
        for bi in block_idx:
            if rng.random() < 0.5:
                flip[bi] *= -1.0
        cnt = int(((xs * flip) * ys > 0).sum())
        if cnt >= obs:
            ge += 1
        flip = np.ones(n)
    p = (ge + 1) / (B + 1)
    return obs, n, p, sizes

# --- 基因间平均相关（CAMERA 式 ρ̄，在两全组织队列的 z 矩阵上） ---
_ns = {"__name__": "_mfdefs", "__file__": os.path.join(HERE, "manuscript_facts.py")}
_src = io.open("manuscript_facts.py", encoding="utf-8").read()
exec(_src[:_src.index('print("载入队列")')], _ns)
load_cohort = _ns["load_cohort"]

def rho_bar(genes, gse):
    c = load_cohort(gse)
    zm = c["zmap_adj"]
    avail = [g for g in genes if g in zm]
    M = np.array([zm[g][:len(c["labels"])] for g in avail])
    M = M - M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, keepdims=True)
    ok = sd[:, 0] > 0
    M = M[ok] / sd[ok]
    R = (M @ M.T) / M.shape[1]
    iu = np.triu_indices(len(avail), 1)
    return float(R[iu].mean()), len(avail)

out = {"UP": {}, "DN": {}}
for key, genes, xs, ys in (("UP", UP_g, UP_x, UP_y), ("DN", DN_g, DN_x, DN_y)):
    obs, n, p_blk, sizes = block_flip_pvalue(genes, xs, ys)
    print("\n== %s：观测 %d/%d" % (key, obs, n))
    print("   块结构：最大块 %d 基因，共 %d 块（top5: %s）" % (sizes[0], len(sizes), sizes[:5]))
    print("   模块分块符号翻转（B=20000）P = %.4f" % p_blk, flush=True)
    rhos = {}
    for gse in ("GSE4183", "GSE117606"):
        r_, nav = rho_bar(genes, gse)
        rhos[gse] = {"rho_bar": round(r_, 4), "n_genes": nav, "n_samples": None}
        print("   ρ̄(%s) = %.3f（%d 基因）" % (gse, r_, nav), flush=True)
    r_avg = (rhos["GSE4183"]["rho_bar"] + rhos["GSE117606"]["rho_bar"]) / 2
    vif = 1 + (n - 1) * r_avg
    z0 = (obs / n - 0.5) / math.sqrt(0.25 / n)
    z_adj = z0 / math.sqrt(vif)
    p_vif = 2 * (1 - 0.5 * (1 + math.erf(abs(z_adj) / math.sqrt(2))))
    print("   VIF = 1+(n-1)ρ̄ = %.2f  | z %.2f → %.2f | 校正二项 P = %.4f"
          % (vif, z0, z_adj, p_vif), flush=True)
    out[key] = {"obs": obs, "n": n, "p_block_flip": round(p_blk, 4),
                "n_blocks": len(sizes), "max_block": sizes[0],
                "rho_bar": {g: rhos[g]["rho_bar"] for g in rhos},
                "vif": round(vif, 2), "z_raw": round(z0, 2),
                "z_vif_adj": round(z_adj, 2), "p_vif_binomial": round(p_vif, 4)}

# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
verdict = all(v["p_block_flip"] < 0.05 and v["p_vif_binomial"] < 0.05 for v in out.values())
print("结论：两法校正后 P %s 0.05 → 方向一致率%s共表达校正后仍显著"
      % ("<" if verdict else "不全部 <", " 均" if verdict else " "))
json.dump({"meta": {"method": "① WGCNA-module block sign-flip (B=20000, blocks flip jointly, preserving co-expression); ② CAMERA-style variance inflation 1+(n-1)*mean_r, binomial z adjusted"},
           "results": out},
          io.open("_cam_gene_correlation_test.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("→ _cam_gene_correlation_test.json")
