# -*- coding: utf-8 -*-
"""补充分析 ⑧：GSE8671 独立验证（B 类遗留 M1）

审稿人诉求（模拟审稿 M1）：在独立数据集上复现方向。
GSE8671（Sabates-Bellver 2007, Gut；PMID 18171984）：**32 例配对**腺瘤 vs 正常黏膜
（无癌组）——正好独立检验「正常 → 腺瘤」这一段的三层模式起点：
  上调分量（间质/髓系激活）在腺瘤期已升高（AUC(腺瘤 vs 正常) > 0.5）
  下调分量（上皮分化/功能）在腺瘤期已下降（AUC < 0.5）
口径与稿件完全一致（无向 mean-z、逐基因队内 z、全局漂移校正 zmap_adj）。
用法：python _gse8671_validation.py   （需联网拉 series matrix，约 1-3 分钟；GPL570 注释走本地缓存）
输出：_gse8671_validation.json
"""
import csv as _csv
import io, json, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

# ---- 复用 manuscript_facts 的通用 load_cohort（内部调 sig_validate.fetch/parse_series） ----
_src = io.open("manuscript_facts.py", encoding="utf-8").read()
_ns = {"__name__": "_mfdefs", "__file__": os.path.join(HERE, "manuscript_facts.py")}
exec(_src[:_src.index('print("载入队列")')], _ns)
load_cohort = _ns["load_cohort"]
PANELS = _ns["PANELS"]
# GSE8671 组织字符串比老表宽（normal mucosa / non-neoplastic …），扩映射防 "other"
_ns["T2K"] = {"tumor": "carcinoma", "adenoma": "adenoma", "ssa": "ssa",
              "adjacent mucosa": "normal", "normal mucosa": "normal",
              "normal colon": "normal", "normal": "normal",
              "normal colonic mucosa": "normal", "normal colon mucosa": "normal",
              "non-neoplastic colon": "normal", "non-neoplastic mucosa": "normal"}

# ---- ★ GSE8671 专用：从 family.soft.gz 载入（matrix 卷的 series matrix 服务端损坏） ----
# 82,197 B 后 gzip 流坏死（FTP/HTTP 双协议、多偏移一致复现 → 服务端文件损坏）。
# soft/ 卷健康 → SOFT 文件含每 GSM 完整表达表，解析口径与 parse_series 对齐。
# ★ manuscript_facts exec sig_validate 时把 defs 存到嵌套的 _ns["_ns"] 里
_sv = _ns["_ns"]
classify3 = _sv["classify3"]
import gzip as _gzip
SOFTF = os.path.join("cohort_hunt", "geo_cache", "GSE8671_family.soft.gz")


def _parse_soft_gse8671(gse):
    samples, cur, in_tbl = [], None, False
    val_col = 1
    with _gzip.open(SOFTF, "rt", encoding="utf-8", errors="replace") as g:
        for ln in g:
            ln = ln.rstrip("\n").rstrip("\r")
            if ln.startswith("^SAMPLE"):
                if cur is not None:
                    samples.append(cur)
                cur = {"title": "", "gpl": None, "chars": [], "expr": {}}
                in_tbl = False
                continue
            if cur is None:
                continue
            if ln.startswith("!Sample_title"):
                cur["title"] = ln.split("=", 1)[1].strip().strip('"')
            elif ln.startswith("!Sample_platform_id"):
                cur["gpl"] = ln.split("=", 1)[1].strip().strip('"')
            elif ln.startswith("!Sample_characteristics_ch1"):
                cur["chars"].append(ln.split("=", 1)[1].strip().strip('"'))
            elif ln.startswith("!sample_table_begin"):
                in_tbl = True
            elif ln.startswith("!sample_table_end"):
                in_tbl = False
            elif in_tbl and ln:
                if ln.startswith("#"):
                    if ln.startswith("#VALUE"):
                        # VALUE 列位置由 "#VALUE = ..." 行序决定，默认第 1 数据列
                        pass
                    continue
                parts = ln.split("\t")
                pid = parts[0].strip().strip('"')
                try:
                    cur["expr"][pid] = float(parts[val_col])
                except (ValueError, IndexError):
                    cur["expr"][pid] = None
    if cur is not None:
        samples.append(cur)

    titles = [s["title"] for s in samples]
    gpl = samples[0]["gpl"]
    labels = []
    for s in samples:
        lab = "other"
        for ch in s["chars"]:
            if "tissue:" in ch.lower():
                v = ch.split(":", 1)[1].strip().lower()
                lab = _ns["T2K"].get(v, "other")
                break
        if lab == "other":
            lab = classify3(s["title"], "", " ".join(s["chars"]))
        labels.append(lab)
    pids = sorted({p for s in samples for p in s["expr"]})
    expr = {p: [s["expr"].get(p) for s in samples] for p in pids}
    print("SOFT: %d 样本, 平台 %s, 探针 %d, 标签 %s"
          % (len(samples), gpl, len(pids),
             {k: labels.count(k) for k in set(labels)}), flush=True)
    print("原始 characteristics:", flush=True)
    for s in samples[:3]:
        print("   ", s["title"], "|", s["chars"], flush=True)
    return {"labels": labels, "expr": expr, "gpl": gpl, "titles": titles}


_ns["parse_series"] = _parse_soft_gse8671
# labels_tissue 会去读损坏的 matrix 缓存 → 屏蔽（labels 已在 SOFT 解析时定好）
_ns["labels_tissue"] = lambda gse: None

SIGINFO = {r["gene"].upper(): r for r in
           _csv.DictReader(io.open("consensus_signature_tier.csv", encoding="utf-8-sig"))}
UP = [g for g, r in SIGINFO.items() if float(r["mean_logFC"]) > 0]
DN = [g for g, r in SIGINFO.items() if float(r["mean_logFC"]) < 0]
ALL313 = sorted(SIGINFO)

print("拉取并载入 GSE8671（首次需联网）…", flush=True)
c = load_cohort("GSE8671")
assert c is not None, "GSE8671 载入失败"
cnt = Counter(c["labels"])
print("  平台 %s | 基因 %d | 标签分布 %s" % (c["gpl"], len(c["zmap"]), dict(cnt)), flush=True)

lab = c["labels"]
iA = [i for i, l in enumerate(lab) if l == "adenoma"]
iN = [i for i, l in enumerate(lab) if l == "normal"]
assert len(iA) >= 20 and len(iN) >= 20, "标签分布异常: %s" % dict(cnt)

import numpy as np


def auc_fast(pos, neg):
    comb = np.concatenate([pos, neg])
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
    n = len(pos)
    return (ranks[:n].sum() - n * (n + 1) / 2.0) / (n * len(neg))


zmap = c["zmap_adj"]
SETSETS = {
    "signature_313": ALL313,
    "UP_component": UP,
    "DN_component": DN,
    "epi_diff_core8": PANELS["上皮分化·核心8"],
    "CAF_panel": PANELS["间质·激活/CAF"],
    "pure_epi_baseline7": PANELS["纯上皮基线·7"],
}
res = {"meta": {"gse": "GSE8671", "gpl": c["gpl"], "n_genes": len(c["zmap"]),
                "labels": dict(cnt), "design": "32 paired adenoma vs normal mucosa "
                "(Sabates-Bellver 2007, PMID 18171984); no carcinoma group",
                "caliber": "undirected mean-z, drift-corrected (zmap_adj); AUC = P(adenoma > normal)"},
       "auc_adenoma_vs_normal": {}}
print("\n无向 AUC（腺瘤 vs 正常；>0.5 = 腺瘤期更高）：")
for key, genes in SETSETS.items():
    gs = [g for g in genes if g in zmap]
    s = {i: np.mean([zmap[g][i] for g in gs]) for i in range(len(lab))}
    a = np.array([s[i] for i in iA])
    n = np.array([s[i] for i in iN])
    v = auc_fast(a, n)
    res["auc_adenoma_vs_normal"][key] = {"n_genes": len(gs), "auc": round(float(v), 4)}
    print("  %-20s n基因=%3d  AUC = %.3f" % (key, len(gs), v), flush=True)

# 逐基因方向复现率（对签名预期方向）
deltas = {}
for g in ALL313:
    z = zmap.get(g)
    if z is None:
        continue
    d = float(np.mean([z[i] for i in iA]) - np.mean([z[i] for i in iN]))
    deltas[g] = d
up_d = [deltas[g] for g in UP if g in deltas]
dn_d = [deltas[g] for g in DN if g in deltas]
gw = [float(np.mean([z[i] for i in iA]) - np.mean([z[i] for i in iN]))
      for z in (zmap[g] for g in zmap)]
p0 = float(np.mean(np.array(gw) > 0))
res["gene_direction"] = {
    "UP_expected_up_frac": round(float(np.mean(np.array(up_d) > 0)), 4),
    "DN_expected_down_frac": round(float(np.mean(np.array(dn_d) < 0)), 4),
    "genome_wise_up_frac_p0": round(p0, 4),
    "n_up_tested": len(up_d), "n_dn_tested": len(dn_d),
}
print("\n逐基因方向复现：UP 分量在腺瘤期上调 %.1f%%（全基因组基线 %.1f%%）| "
      "DN 分量在腺瘤期下调 %.1f%%"
      % (100 * res["gene_direction"]["UP_expected_up_frac"],
         100 * p0, 100 * res["gene_direction"]["DN_expected_down_frac"]))

# ---- 补充统计：二项检验（对全基因组基线）+ delta 与 mean_logFC 的 Spearman ----
from math import erf, sqrt

def binom_cdf(k, n, p):
    # 精确二项 CDF（小 n 直接递推）
    from math import comb
    return sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(0, k + 1))

k_up = int(round(res["gene_direction"]["UP_expected_up_frac"] * len(up_d)))
n_up = len(up_d)
p_up_below = binom_cdf(k_up, n_up, p0)          # UP 基因「上调率显著低于基线」
k_dn = int(round(res["gene_direction"]["DN_expected_down_frac"] * len(dn_d)))
n_dn = len(dn_d)
p0_dn = 1 - p0
p_dn_above = 1 - binom_cdf(k_dn - 1, n_dn, p0_dn)  # DN 下调率「显著高于随机下调基线」
res["gene_direction"]["UP_vs_baseline_binom_P"] = float("%.3g" % max(p_up_below, 1e-300))
res["gene_direction"]["DN_vs_baseline_binom_P"] = float("%.3g" % max(p_dn_above, 1e-300))
print("二项检验：UP 上调率低于基线 P = %.3g | DN 下调率高于随机基线 P = %.3g"
      % (p_up_below, p_dn_above), flush=True)

# Spearman（手写平均秩）
sig_tier = {g: float(SIGINFO[g]["mean_logFC"]) for g in deltas if g in SIGINFO}
xs = np.array([deltas[g] for g in sig_tier])
ys = np.array([sig_tier[g] for g in sig_tier])

def rankdata(a):
    order = a.argsort(kind="mergesort")
    r = np.empty(len(a))
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0 + 1
        i = j + 1
    return r

rx, ry = rankdata(xs), rankdata(ys)
rho = float(np.corrcoef(rx, ry)[0, 1])
res["gene_direction"]["spearman_delta_vs_meanlogFC"] = round(rho, 4)
res["gene_direction"]["n_gene_corr"] = len(sig_tier)
print("Spearman（GSE8671 Δ 与签名 mean_logFC）r = %.3f (n=%d)" % (rho, len(sig_tier)),
      flush=True)

# ---- 对照：内部三个队列同口径的 adenoma vs normal（结论对照表用） ----
try:
    inhouse = json.load(io.open("_adenoma_vs_normal_inhouse.json", encoding="utf-8"))
    res["inhouse_reference"] = inhouse
except Exception:
    pass

json.dump(res, io.open("_gse8671_validation.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\n→ _gse8671_validation.json")
