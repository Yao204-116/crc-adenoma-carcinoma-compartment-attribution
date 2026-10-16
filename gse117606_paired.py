# -*- coding: utf-8 -*-
"""
GSE117606 患者配对验证（第三独立队列 / 第二平台）
================================================================================
为什么为它单独写一个脚本：

  GSE117606 = 208 例，设计是「同一患者的 A 腺瘤 / B 癌 / C 癌旁黏膜」：
      Tumor 74 / Adjacent Mucosa 65 / Adenoma 59 / SSA 10
      患者组合：Adenoma+ACM+Tumor 48 人，Adenoma+Tumor 8 人 → **56 组腺瘤-癌配对**
  配对设计能消掉个体间变异（年龄/性别/遗传背景/菌群/饮食/取材批次），
  这是检验「组成效应 vs 细胞内在」时最有力的设计——比多队列更有说服力。

  另外两件事必须在本脚本里纠正：
   ① **剔除 SSA**（无蒂锯齿状腺瘤）。上一轮分类器把 10 例 SSA 算进了腺瘤
      （所以日志里印的是「腺瘤 69」）。SSA 走锯齿状通路（BRAF/CIMP），分子特征
      更接近癌，混进腺瘤会**低估**签名区分力，且与我们既定的「排除锯齿状通路」冲突。
   ② 该队列 title 里有 3 例 code=A 但描述=Tumor，所以**以 tissue 字段为准**判类型。

输出三块：
  A. 全体非配对（59 腺瘤 vs 74 癌）—— 与前两队列口径一致，可比
  B. ★ 56 组配对：Δ = 癌 − 腺瘤，符号一致率 + Wilcoxon 符号秩检验（手写，无 scipy）
  C. 48 组三元组的「腺瘤 vs 癌旁正常」配对对照
用法：python gse117606_paired.py
"""
import csv, gzip, importlib.util, json, math, os, re, shutil, statistics, sys
from collections import Counter, defaultdict

BASE = r"C:/Users/14003/WorkBuddy/2026-10-08-09-31-50/geo_crc/validation_search"
CACHE = os.path.join(BASE, "cohort_hunt", "geo_cache")
RES_JSON = os.path.join(BASE, "cohort_hunt", "sig_validation_results.json")
BAK_JSON = os.path.join(BASE, "sig_validation_results_4cohort.json")
OUT_JSON = os.path.join(BASE, "gse117606_paired_results.json")

# ---------------------------------------------------------------------------
# 加载 sig_validate 的既有函数（会复用 parse_series / load_annot / load_entrez_map / auc / fmt_auc）
# ★ 坑：sig_validate.py 是脚本式写法，模块级末尾有 `results={}` + `json.dump(...)`，
#   直接 import 会把已经跑好的 4 队列结果覆盖成 `{}`。所以这里先备份、后恢复。
#   同时把 sys.argv 清空，避免它模块级的 `for gse in sys.argv[1:]` 主循环被误触发。
# ---------------------------------------------------------------------------
if os.path.exists(RES_JSON):
    shutil.copy2(RES_JSON, BAK_JSON)

sys.path.insert(0, BASE)
_argv = sys.argv[:]
sys.argv = [sys.argv[0]]
try:
    _spec = importlib.util.spec_from_file_location("_sv", os.path.join(BASE, "sig_validate.py"))
    sv = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(sv)
finally:
    sys.argv = _argv
    if os.path.exists(BAK_JSON):
        shutil.copy2(BAK_JSON, RES_JSON)   # 恢复被覆盖的 4 队列结果

UP, DN, SIG = sv.UP, sv.DN, sv.SIG
PANELS = sv.PANELS

# 追加面板：CLCA1 不在原 PANELS 里，而它恰恰是 LCM 里掉得最狠的那个基因（16.71→5.59）
EXTRA = {
    "上皮分化·核心8": "CLCA1 ZG16 ITLN1 MUC2 REG4 FCGBP SPINK4 SI",
    "纯上皮基线·7": "EPCAM KRT8 KRT18 KRT19 CDH1 CDH17 VIL1",
    "泛间质·6": "FAP THBS2 COL11A1 INHBA POSTN ASPN",
    "髓系·激活": ("SPP1 CCL18 CCL2 CCL3 CCL4 CCL20 APOE TREM2 MMP9 MMP12 SIGLEC1 CXCL8 "
                  "IL10 IDO1 CD274 PDCD1LG2 LGALS3 S100A8 S100A9 S100A12 FCN1 C5AR1 "
                  "MNDA FCGR3B MMP8 TCN1 CCL19"),
}
EXTRA = {k: set(v.split()) for k, v in EXTRA.items()}
ALLPANELS = dict(PANELS)
ALLPANELS.update(EXTRA)

TISSUE2KIND = {
    "tumor": "carcinoma",
    "adenoma": "adenoma",
    "ssa": "ssa",                 # ★ 锯齿状通路，剔除
    "adjacent mucosa": "normal",
}


# ===========================================================================
# 手写统计（环境无 numpy/scipy）
# ===========================================================================
def wilcoxon_signed(deltas):
    """
    Wilcoxon 符号秩检验（剔除 Δ=0 的对）。返回 (W+, W-, n, z, p_two_sided)。
    n 较大时用正态近似（含 ties 的平均秩修正），与本场景（n≈56）匹配。
    """
    dz = [x for x in deltas if x != 0]
    n = len(dz)
    if n == 0:
        return None
    order = sorted(range(n), key=lambda k: abs(dz[k]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(dz[order[j + 1]]) == abs(dz[order[i]]):
            j += 1
        r = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    Wp = sum(ranks[k] for k in range(n) if dz[k] > 0)
    Wm = sum(ranks[k] for k in range(n) if dz[k] < 0)
    mu = n * (n + 1) / 4.0
    sd = math.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    z = (Wp - mu - 0.5) / sd if sd else 0.0     # 连续性校正 0.5
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return Wp, Wm, n, z, p


def sign_test(n_pos, n_tot):
    """符号检验（精确二项，双侧）。用于「有多少对是同向的」。"""
    if n_tot == 0:
        return 1.0
    k = max(n_pos, n_tot - n_pos)
    p = sum(math.comb(n_tot, i) for i in range(k, n_tot + 1)) / (2.0 ** n_tot)
    return min(1.0, 2 * p)


def mw_p(pos, neg):
    """
    Mann-Whitney U 的双侧渐近 p（正态近似 + 连续性校正）。
    n<8 时渐近近似偏乐观，小队列的 p 只当参考——所以报告里小队列标 n 而不强调 p。
    """
    n1, n2 = len(pos), len(neg)
    if n1 < 2 or n2 < 2:
        return None
    U = sv.auc(pos, neg) * n1 * n2
    mu = n1 * n2 / 2.0
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = (U - mu - 0.5) / sd if sd else 0.0
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def fnum(x, w=3):
    return ("%+." + str(w) + "f") % x if isinstance(x, (int, float)) else "—"


def fpa(p):
    if p is None:
        return "—"
    if p < 1e-4:
        return "%.2e" % p
    return "%.4f" % p


# ===========================================================================
# 1. 读 series matrix 头部拿 title / tissue 字段（authoritative），并解析患者配对
# ===========================================================================
def head_tissue(gse):
    p = os.path.join(CACHE, gse + "_series_matrix.txt.gz")
    titles, tissue, gsms = [], [], []
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!Sample_geo_accession"):
                gsms = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!Sample_characteristics_ch1") and "tissue:" in ln:
                tissue = [x.strip().strip('"').replace("tissue:", "").strip()
                          for x in ln.rstrip("\n").split("\t")[1:]]
            elif ln.startswith("!series_matrix_table_begin"):
                break
    return titles, tissue, gsms


print("=" * 100)
print("GSE117606 患者配对验证")
print("=" * 100)
titles, tissue, gsms = head_tissue("GSE117606")
print("  series matrix 头部：title %d / tissue %d / gsm %d" % (len(titles), len(tissue), len(gsms)))

recs = []
for i, (t, ti) in enumerate(zip(titles, tissue)):
    m = re.match(r"^(\S+?):\s*(.+?)\s+of\s+Patient\s+(\S+)$", t)
    pid = m.group(3) if m else "?"
    kind = TISSUE2KIND.get(ti.lower(), "other")
    recs.append({"i": i, "pid": pid, "kind": kind, "tissue": ti,
                 "gsm": gsms[i] if i < len(gsms) else ""})

kc = Counter(r["kind"] for r in recs)
print("\n  样本清点（按 tissue 字段，权威）：")
print("    肿瘤 %d | 腺瘤 %d | 癌旁正常 %d | ★SSA(锯齿状，剔除) %d"
      % (kc["carcinoma"], kc["adenoma"], kc["normal"], kc["ssa"]))
print("    （对比上一轮分类器：腺瘤报 69 = 普通腺瘤 59 + SSA 10 → 本轮已拆开）")

byp = defaultdict(lambda: defaultdict(list))
for r in recs:
    byp[r["pid"]][r["kind"]].append(r["i"])

pairs_ac = sorted(p for p, d in byp.items() if d["adenoma"] and d["carcinoma"])
pairs_all3 = sorted(p for p, d in byp.items()
                    if d["adenoma"] and d["carcinoma"] and d["normal"])
print("\n  配对结构：")
print("    腺瘤-癌配对人 %d（其中同时有癌旁正常的三元组 %d）" % (len(pairs_ac), len(pairs_all3)))

# ===========================================================================
# 2. 注释 + 表达矩阵 + z（逻辑与 sig_validate 主循环一致，保证与前两队列同口径）
# ===========================================================================
d = sv.parse_series("GSE117606")
if d is None:
    sys.exit("!! 无法解析 series matrix")
ann = sv.load_annot(d["gpl"]) or sv.load_entrez_map(d["expr"].keys())
print("  平台 %s | 注释探针 %d 条" % (d["gpl"], len(ann)))

gene_rows = {}
for pid_, vals in d["expr"].items():
    g = ann.get(pid_)
    if not g:
        continue
    good = [v for v in vals if v is not None]
    if len(good) < len(vals) * 0.8:
        continue
    m = statistics.mean(good)
    if g not in gene_rows or m > gene_rows[g][0]:
        gene_rows[g] = (m, vals)
print("  可比对基因 %d 个" % len(gene_rows))

NG = len(d["labels"])
sample_med = []
for i in range(NG):
    v = [gene_rows[g][1][i] for g in gene_rows
         if i < len(gene_rows[g][1]) and gene_rows[g][1][i] is not None]
    if v:
        sample_med.append(statistics.median(v))
if sample_med and statistics.median(sample_med) > 30:
    print("  检测到线性尺度 → 做 log2(x+1)")
    for g in gene_rows:
        gene_rows[g] = (gene_rows[g][0], [math.log2(v + 1) if v is not None else None
                                          for v in gene_rows[g][1]])


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

# 基因名 → 可用性检查（关键基因是否被平台测到）
KEYCHK = ["CLCA1", "FAP", "THBS2", "EPCAM", "KRT18", "MUC2", "SPP1", "CDH1"]
print("\n  关键基因是否可比对：%s"
      % "  ".join("%s=%s" % (g, "有" if g in zmap else "★缺") for g in KEYCHK))


def score(i, genes):
    zz = [zmap[g][i] for g in genes if g in zmap and i < len(zmap[g])]
    return statistics.mean(zz) if zz else None


def score_sig(i):
    """
    ★ 全签名的正确定义 = 上调分量均值 − 下调分量均值（与 sig_validate.py 的 score() 一致）。
    踩坑：本脚本初版误写成「313 个基因的 z 均值」（score(i, UP+DN)），
      由于上下调基因互相抵消，得到的是完全不同的量（癌中位 +0.114 vs 正确 +1.194）。
      判据：逐基因方向符合率仍报 99.0%（那条逻辑是对的），但全签名 AUC 对不上 → 定位到此。
    """
    u, d_ = score(i, UP), score(i, DN)
    if u is None or d_ is None:
        return None
    return u - d_


# 统一的「签名/分量」定义，A/B/C/D 四段共用，避免各段口径漂移
SIGSETS = [
    ("全签名", score_sig, UP + DN),
    ("上调分量(间质/髓系激活)", lambda i: score(i, UP), UP),
    ("下调分量(上皮功能关闭)", lambda i: score(i, DN), DN),
]


def auc_wrap(pos, neg):
    return sv.auc(pos, neg)


def fmt(a, p="癌", n="腺瘤"):
    return sv.fmt_auc(a, p, n)


# ===========================================================================
# 3. 分析 A：全体非配对（剔除 SSA 后）
# ===========================================================================
idx_ad = [r["i"] for r in recs if r["kind"] == "adenoma"]
idx_ca = [r["i"] for r in recs if r["kind"] == "carcinoma"]
idx_no = [r["i"] for r in recs if r["kind"] == "normal"]

print("\n" + "-" * 100)
print("A. 全体非配对（剔除 10 例 SSA 后）：腺瘤 n=%d vs 癌 n=%d" % (len(idx_ad), len(idx_ca)))
print("-" * 100)

resA = {}
for tag, fn, genes in SIGSETS:
    a = [fn(i) for i in idx_ad]
    c = [fn(i) for i in idx_ca]
    a = [x for x in a if x is not None]
    c = [x for x in c if x is not None]
    aa = auc_wrap(c, a)
    pp = mw_p(c, a)
    resA[tag] = {"auc": aa, "p": pp}
    print("  %-24s n基因=%3d  AUC = %-22s p=%-9s 腺瘤中位 %s / 癌中位 %s"
          % (tag, len([g for g in genes if g in zmap]), fmt(aa), fpa(pp), fnum(statistics.median(a)),
             fnum(statistics.median(c))))
print("  逐基因方向符合率：", end="")
hit = tot = 0
for r in SIG:
    g = r["gene"]
    if g not in zmap:
        continue
    tot += 1
    exp = 1 if float(r["mean_logFC"]) > 0 else -1
    obs = 1 if statistics.mean([zmap[g][i] for i in idx_ca]) > \
               statistics.mean([zmap[g][i] for i in idx_ad]) else -1
    hit += (exp == obs)
print("%d / %d = %.1f%%" % (hit, tot, 100.0 * hit / max(tot, 1)))

# A2. 小面板单独 AUC —— 面板越小越容易落地成抗体，所以必须单看
print("\n  ★ 关键面板的非配对 AUC（面板越小越值得做成 mIF 抗体组合）")
resA2 = {}
for k in ["上皮分化·核心8", "泛间质·6", "纯上皮基线·7", "间质·激活/CAF", "间质·基线"]:
    genes = ALLPANELS[k]
    a = [score(i, genes) for i in idx_ad]
    c = [score(i, genes) for i in idx_ca]
    a = [x for x in a if x is not None]
    c = [x for x in c if x is not None]
    aa = auc_wrap(c, a)
    resA2[k] = {"auc": aa, "p": mw_p(c, a), "n": len([g for g in genes if g in zmap])}
    print("     %-16s n=%2d  AUC = %-22s p=%s"
          % (k, resA2[k]["n"], fmt(aa), fpa(resA2[k]["p"])))

# ===========================================================================
# 4. 分析 B：★ 56 组患者配对（消掉个体间变异）
# ===========================================================================
print("\n" + "-" * 100)
print("B. ★ 患者配对：%d 组（Δ = 癌 − 腺瘤，单位=队列内 z；符号应全部为正）" % len(pairs_ac))
print("-" * 100)

resB = {}
for tag, fn, genes in SIGSETS:
    deltas, ndet = [], 0
    for p in pairs_ac:
        ia = byp[p]["adenoma"][0]
        ic = byp[p]["carcinoma"][0]
        sa, sc = fn(ia), fn(ic)
        if sa is None or sc is None:
            ndet += 1
            continue
        deltas.append(sc - sa)
    npos = sum(1 for x in deltas if x > 0)
    wt = wilcoxon_signed(deltas)
    st = sign_test(npos, len(deltas))
    resB[tag] = {"n_pair": len(deltas), "n_pos": npos,
                 "median_delta": statistics.median(deltas),
                 "wilcoxon_p": wt[4] if wt else None,
                 "sign_test_p": st}
    print("  %-24s 配对 %d | 癌>腺瘤 %d/%d = %5.1f%% | 中位 Δ %s | Wilcoxon p=%s | 符号检验 p=%s"
          % (tag, len(deltas), npos, len(deltas), 100.0 * npos / max(len(deltas), 1),
             fnum(statistics.median(deltas)), fpa(wt[4] if wt else None), fpa(st)))

# 面板级配对
print("\n  ★ 面板级配对（同一患者内 癌 vs 腺瘤）")
resB_panel = {}
for k, genes in ALLPANELS.items():
    deltas = []
    for p in pairs_ac:
        sa = score(byp[p]["adenoma"][0], genes)
        sc = score(byp[p]["carcinoma"][0], genes)
        if sa is None or sc is None:
            continue
        deltas.append(sc - sa)
    if len(deltas) < 5:
        continue
    npos = sum(1 for x in deltas if x > 0)
    wt = wilcoxon_signed(deltas)
    resB_panel[k] = {"n": len(deltas), "n_pos": npos, "median_delta": statistics.median(deltas),
                     "wilcoxon_p": wt[4] if wt else None}
    arrow = "↑癌" if npos > len(deltas) / 2 else "↓癌"
    print("     %-16s n基因=%3d | 配对 %2d | 癌>腺瘤 %2d/%2d = %5.1f%% %s | 中位 Δ %s | p=%s"
          % (k, len([g for g in genes if g in zmap]), len(deltas), npos, len(deltas),
             100.0 * npos / len(deltas), arrow,
             fnum(statistics.median(deltas)), fpa(wt[4] if wt else None)))

# ===========================================================================
# 5. 分析 C：48 组三元组的「腺瘤 vs 癌旁正常」配对
# ===========================================================================
print("\n" + "-" * 100)
print("C. 配对对照：腺瘤 vs 癌旁正常（%d 组三元组）—— 判断签名在序贯的哪一步开始变"
      % len(pairs_all3))
print("-" * 100)

resC = {}
for tag, fn, genes in SIGSETS:
    d1 = []   # 腺瘤 − 正常
    for p in pairs_all3:
        sn = fn(byp[p]["normal"][0])
        sa = fn(byp[p]["adenoma"][0])
        if sn is None or sa is None:
            continue
        d1.append(sa - sn)
    npos = sum(1 for x in d1 if x > 0)
    wt = wilcoxon_signed(d1)
    resC[tag] = {"n": len(d1), "n_pos": npos, "median_delta": statistics.median(d1) if d1 else None,
                 "wilcoxon_p": wt[4] if wt else None}
    print("  %-24s 配对 %d | 腺瘤>正常 %d/%d = %5.1f%% | 中位 Δ(腺瘤−正常) %s | Wilcoxon p=%s"
          % (tag, len(d1), npos, len(d1), 100.0 * npos / max(len(d1), 1),
             fnum(statistics.median(d1)) if d1 else "—", fpa(wt[4] if wt else None)))

# 面板级
print("\n  ★ 面板级配对（腺瘤 − 癌旁正常）")
resC_panel = {}
for k, genes in ALLPANELS.items():
    d1 = []
    for p in pairs_all3:
        sn = score(byp[p]["normal"][0], genes)
        sa = score(byp[p]["adenoma"][0], genes)
        if sn is None or sa is None:
            continue
        d1.append(sa - sn)
    if len(d1) < 5:
        continue
    npos = sum(1 for x in d1 if x > 0)
    wt = wilcoxon_signed(d1)
    resC_panel[k] = {"n": len(d1), "n_pos": npos, "median_delta": statistics.median(d1),
                     "wilcoxon_p": wt[4] if wt else None}
    print("     %-16s 配对 %2d | 腺瘤>正常 %2d/%2d = %5.1f%% | 中位 Δ %s | p=%s"
          % (k, len(d1), npos, len(d1), 100.0 * npos / len(d1),
             fnum(statistics.median(d1)), fpa(wt[4] if wt else None)))

# ===========================================================================
# 5b. ★ 必须解释的反常：剔除 SSA 后非配对「全签名」AUC 从 0.918 掉到 0.830，
#     而两个分量 AUC 几乎不变（0.901→0.906、0.910→0.913）。SSA 到底落在哪？
# ===========================================================================
idx_ssa = [r["i"] for r in recs if r["kind"] == "ssa"]
print("\n" + "-" * 100)
print("D. SSA 组（n=%d，锯齿状通路）在签名分上的位置 —— 解释「剔 SSA 后 AUC 为何下降」"
      % len(idx_ssa))
print("-" * 100)
resD = {}
for tag, fn, genes in SIGSETS:
    row = {}
    for nm, idxs in [("正常", idx_no), ("腺瘤", idx_ad), ("★SSA", idx_ssa), ("癌", idx_ca)]:
        vals = [fn(i) for i in idxs]
        vals = [x for x in vals if x is not None]
        row[nm] = statistics.median(vals) if vals else None
    resD[tag] = row
    print("  %-24s 正常 %s | 腺瘤 %s | ★SSA %s | 癌 %s"
          % (tag, fnum(row["正常"]), fnum(row["腺瘤"]), fnum(row["★SSA"]), fnum(row["癌"])))
print("  读法：若 SSA 中位数与「正常」同侧、离「癌」很远，则它是把腺瘤组整体拉向低分的一群，")
print("        剔除后腺瘤组分数上移 → 与癌重叠变多 → 全签名 AUC 下降。分量 AUC 稳，说明分量更稳健。")

# ===========================================================================
# 6. 落盘
# ===========================================================================
out = {
    "gse": "GSE117606", "gpl": d["gpl"],
    "counts": dict(kc), "n_pairs_ac": len(pairs_ac), "n_pairs_all3": len(pairs_all3),
    "A_unpaired": resA, "A_panel_auc": resA2,
    "A_gene_concordance": hit / max(tot, 1), "A_n_tested": tot,
    "B_paired": resB, "B_panel_paired": resB_panel,
    "C_adenoma_vs_normal_paired": resC, "C_panel_paired": resC_panel,
    "D_ssa_position": resD,
    "key_genes_present": {g: (g in zmap) for g in KEYCHK},
}
json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), ensure_ascii=False, indent=2, default=str)
print("\n[saved] %s" % OUT_JSON)
print("DONE")
