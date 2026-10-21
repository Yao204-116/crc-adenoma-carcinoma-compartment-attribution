# -*- coding: utf-8 -*-
"""
三层对比（正常 / 腺瘤 / 癌）在全部可用队列上的双层签名表现。

为什么要新增这个口径：
  前几轮所有 AUC 都是「腺瘤 vs 癌」——但用户的实验室**没有腺瘤组织**，
  唯一可能落地的 mIF 验证是「癌 vs 癌旁正常」。
  所以必须知道：双层签名在「癌 vs 正常」上还成不成立。
  如果不成立，mIF 这条路连尾巴都留不下。

★ 同时检查「癌旁正常」是否已被场效应污染（癌旁 与 腺瘤 的相似度）。
  若癌旁比真正的正常更靠近腺瘤，说明「癌旁=正常」这个假设本身有问题 ——
  这会让所有用「癌 vs 癌旁」的 CRC 研究都带上系统偏差。

复用 sig_validate 的解析与注释逻辑（模块级代码会重写结果文件，故仅取函数）。
"""
import json, os, re, sys, csv, math, statistics, importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)

# ── 动态载入 sig_validate（只取它的解析/注释函数） ──
spec = importlib.util.spec_from_file_location("svmod", os.path.join(HERE, "sig_validate.py"))
sv = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(sv)
except SystemExit:
    pass
except Exception as e:
    print("模块载入异常（可能已在 __main__ 里跑完）:", type(e).__name__, e)


def load_cohort(gse):
    """★ 从 collide_all.py 搬过来（它在 sv 里不存在）。
    返回 {"d": 含 labels, "zmap": 基因→z 行, "gpl"}。"""
    d = sv.parse_series(gse)
    if d is None:
        return None
    ann = sv.load_annot(d["gpl"]) or sv.load_entrez_map(d["expr"].keys())
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
    if med and statistics.median(med) > 30:      # 未 log2 → 现场补
        for g in gene_rows:
            gene_rows[g] = (gene_rows[g][0],
                            [math.log2(v + 1) if v is not None else None
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
    return {"d": d, "zmap": zmap, "gpl": d["gpl"], "titles": d["titles"]}


print("UP=%d DN=%d" % (len(sv.UP), len(sv.DN)))
P = sv.PANELS if hasattr(sv, "PANELS") else {}

EXTRA = {
    "泛间质·6": "FAP THBS2 COL11A1 INHBA POSTN ASPN",
    "上皮分化·核心8": "CLCA1 ZG16 ITLN1 MUC2 REG4 FCGBP SPINK4 SI",
    "纯上皮基线·7": "EPCAM KRT8 KRT18 KRT19 CDH1 CDH17 VIL1",
    "相代谢·5": "ADH1C CYP3A5 FMO5 NAT2 UGT2B15",
    "潘氏细胞·4": "DEFA5 DEFA6 OLFM4 SPINK1",
    "黏液杯状·5": "CLCA1 MUC2 FCGBP SPINK4 WFDC2",
}
EXTRA = {k: v.split() for k, v in EXTRA.items()}
ALLP = dict(P)
ALLP.update(EXTRA)

SETS = {
    "313全签名": list(sv.UP) + list(sv.DN),
    "上调分量(210)": list(sv.UP),
    "下调分量(103)": list(sv.DN),
}
SETS.update(ALLP)

COHORTS = ["GSE4183", "GSE117606", "GSE15960"]   # 只跑已缓存的，避免重新下载被 SIGTERM


def label_of(idx, labels):
    return labels[idx]


def mean_score(zmap, idx, genes):
    v = [zmap[g][idx] for g in genes if g in zmap and idx < len(zmap[g])]
    return statistics.mean(v) if v else None


def auc(pos, neg):
    if not pos or not neg:
        return None
    n = 0
    for a in pos:
        for b in neg:
            n += 1 if a > b else (0.5 if a == b else 0)
    return n / (len(pos) * len(neg))


def fmt(a, pos="癌", neg="正常"):
    if a is None:
        return "—"
    d = max(a, 1 - a)
    if abs(a - 0.5) < 1e-9:
        return "%.3f (无方向)" % d
    return "%.3f (%s>%s)" % (d, pos if a > 0.5 else neg, neg if a > 0.5 else pos)


RES = {}
# ★ 记录被跳过的队列。一旦非空，本次产物就是**不完整**的，不能覆盖归档文件。
#   旧实现无论跳过几个队列都照写同名文件，于是"跑一遍"会把完整结果换成残缺版本
#   （实测：GSE4183 的 carcinoma 整组消失、GSE117606 整块变 None），而退出码仍是 0。
#   读者拿它对照论文数字，只会得出"跑不出来"。
SKIPPED = []
for gse in COHORTS:
    try:
        c = load_cohort(gse)
    except Exception as e:
        c = None
        print("  %s 载入失败: %s" % (gse, type(e).__name__))
    if not c:
        print("  %s 跳过（无数据/注释）" % gse)
        SKIPPED.append(gse)
        continue
    lab = c["d"]["labels"]          # ★ labels 在 d 里，不在顶层
    zmap = c["zmap"]
    avail = set(zmap)
    idx = {t: [i for i, l in enumerate(lab) if l == t] for t in set(lab)}
    print("\n" + "=" * 78)
    print("%s（%s）" % (gse, c.get("gpl", "?")))
    for t in ("normal", "adenoma", "ssa", "carcinoma", "other"):
        if idx.get(t):
            print("   %-11s n=%d" % (t, len(idx[t])))
    RES[gse] = {"gpl": c.get("gpl"), "n": {t: len(v) for t, v in idx.items() if v}}
    if not (idx.get("normal") and idx.get("carcinoma")):
        print("   → 缺 normal 或 carcinoma，跳过 AUC")
        continue
    RES[gse]["auc"] = {}
    for name, genes in SETS.items():
        gs = [g for g in genes if g in avail]
        if len(gs) < 3:
            continue
        no = [mean_score(zmap, i, gs) for i in idx["normal"]]
        no = [x for x in no if x is not None]
        ca = [mean_score(zmap, i, gs) for i in idx["carcinoma"]]
        ca = [x for x in ca if x is not None]
        row = {"n_genes": len(gs), "auc_ca_vs_no": auc(ca, no),
               "median_normal": statistics.median(no), "median_carcinoma": statistics.median(ca)}
        if idx.get("adenoma"):
            ad = [mean_score(zmap, i, gs) for i in idx["adenoma"]]
            ad = [x for x in ad if x is not None]
            row["auc_ad_vs_no"] = auc(ad, no)
            row["auc_ca_vs_ad"] = auc(ca, ad)
            row["median_adenoma"] = statistics.median(ad)
        RES[gse]["auc"][name] = row
    # ── 打印：重点是 癌 vs 正常 和 癌 vs 腺瘤 的对比 ──
    hdr = ["基因集", "n", "癌vs正常", "腺瘤vs正常", "癌vs腺瘤", "正常中位", "腺瘤中位", "癌中位"]
    print("   %-18s %4s %-20s %-20s %-20s %8s %8s %8s" % tuple(hdr))
    for name, r in RES[gse]["auc"].items():
        print("   %-18s %4d %-20s %-20s %-20s %+8.3f %+8.3f %+8.3f"
              % (name, r["n_genes"],
                 fmt(r["auc_ca_vs_no"]), fmt(r.get("auc_ad_vs_no")), fmt(r.get("auc_ca_vs_ad")),
                 r["median_normal"], r.get("median_adenoma") or 0, r["median_carcinoma"]))

if SKIPPED:
    print("\n" + "!" * 78)
    print("!! 因缺数据被跳过：%s" % ", ".join(SKIPPED))
    print("!! 本次结果**不完整**，故不覆盖归档文件，改写到 *.partial.json")
    print("!! 完整复现请先准备 cohort_hunt/geo_cache/（见 REPRODUCIBILITY.md 9.1）")
    print("!" * 78)
SUF = ".partial" if SKIPPED else ""

json.dump(RES, open(os.path.join(HERE, "three_layer_results%s.json" % SUF), "w",
                    encoding="utf-8"), ensure_ascii=False, indent=1)

# ===========================================================================
# 汇总：三层是否**单调下降**，还是「腺瘤期不降反升、癌期才塌陷」
# —— 这是本轮数据里冒出的唯一反直觉现象，必须跨队列钉死
# ===========================================================================
print("\n" + "=" * 78)
print("★ 三层模式汇总（正常 / 腺瘤 / 癌 的中位 z 分数）")
print("   模式判据： D=单调下降(no>ad>ca)  P=腺瘤峰(ad>=no 且 ad>ca)  A=腺瘤已降(ca>=ad)  O=其他")
print("=" * 78)
WATCH = ["313全签名", "上调分量(210)", "下调分量(103)", "上皮·分化/功能",
         "上皮分化·核心8", "相代谢·5", "黏液杯状·5", "潘氏细胞·4",
         "间质·激活/CAF", "泛间质·6", "纯上皮基线·7"]
SUMMARY = {}
for name in WATCH:
    print("\n── %s ──" % name)
    SUMMARY[name] = {}
    for gse in RES:
        r = RES[gse].get("auc", {}).get(name)
        if not r or "median_adenoma" not in r:
            continue
        no, ad, ca = r["median_normal"], r["median_adenoma"], r["median_carcinoma"]
        if no > ad > ca:
            pat = "D 单调降"
        elif ad >= no and ad > ca:
            pat = "P 腺瘤峰"
        elif ca >= ad:
            pat = "A 腺瘤已降"
        else:
            pat = "O 其他"
        SUMMARY[name][gse] = {"normal": no, "adenoma": ad, "carcinoma": ca, "pattern": pat}
        print("   %-10s 正常 %+7.3f → 腺瘤 %+7.3f → 癌 %+7.3f   %s" % (gse, no, ad, ca, pat))
json.dump(SUMMARY, open(os.path.join(HERE, "three_layer_pattern%s.json" % SUF), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)
print("\n→ three_layer_results%s.json + three_layer_pattern%s.json" % (SUF, SUF))
sys.exit(2 if SKIPPED else 0)
