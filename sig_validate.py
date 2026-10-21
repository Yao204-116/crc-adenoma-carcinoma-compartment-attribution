# -*- coding: utf-8 -*-
"""
独立验证：把 4 队列共识得到的 313 基因签名，拿到"没参与构建"的队列上去测
要点：
  ① 必须用未参与构建的队列（否则是循环论证）
  ② 必须自己数样本分组（不信任 GEO 标题/摘要）
  ③ 报告三件事：签名 AUC、逐基因方向符合率、各细胞面板中位 logFC 是否复现
用法：python sig_validate.py GSE32323 GSE4183 ...
"""
import csv, gzip, io, json, math, os, re, statistics, sys, time, urllib.request, urllib.parse
from collections import Counter

UA = {"User-Agent": "Mozilla/5.0 (research meta-analysis; contact-local)"}
BASE = os.path.dirname(os.path.abspath(__file__))   # 仓库自足：相对脚本位置，不写死本机路径
SIGFILE = os.path.join(BASE, "consensus_signature_tier.csv")
CACHE = os.path.join(BASE, "cohort_hunt", "geo_cache")
os.makedirs(CACHE, exist_ok=True)


def fetch(url, timeout=600):
    """
    返回 **原始压缩字节**（调用方自己 gzip.decompress）。
    ★ BUG 修正：旧版命中缓存时返回 gzip.open(c).read()（已解压），而 parse_series 又解压一次，
      导致「已缓存的文件反而解析失败」。现在无论缓存还是网络都返回原始字节。
    ★ 另：NCBI FTP 单连接实测只有 ~100 KB/s，且会间歇性停顿 —— 所以 timeout 放到 600s，
      并优先使用 geo_dl.py 预先并行下好的缓存。
    """
    cache = os.path.join(CACHE, os.path.basename(url))
    if os.path.exists(cache) and os.path.getsize(cache) > 500:
        return open(cache, "rb").read()
    last = None
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
                d = r.read()
            with open(cache, "wb") as f:
                f.write(d)
            return d
        except Exception as e:
            last = e
            print("     [fetch retry %d] %s : %s" % (i + 1, os.path.basename(url), type(e).__name__), flush=True)
            time.sleep(3 * (i + 1))
    print("     [fetch FAIL] %s : %s" % (os.path.basename(url), last), flush=True)
    return None


def matrix_url(gse):
    num = gse[3:]
    stem = num[:-3] if len(num) > 3 else ""
    return ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%snnn/%s/matrix/%s_series_matrix.txt.gz"
            % (stem, gse, gse))


def gpl_url(gpl):
    num = gpl[3:]
    stem = num[:-3] if len(num) > 3 else ""
    return "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL%snnn/%s/annot/%s.annot.gz" % (stem, gpl, gpl)


NORM = ("normal", "healthy", "non-tumor", "nontumor", "non-neoplastic",
        "adjacent", "control", "uninvolved")

# 复用 cohort_hunt2 里已用 4 个已知队列验证过的「优先级+回退」分类器。
# 注意：**不要**再用下面的朴素子串法 —— 它会把 GSE8671 的正常黏膜判成腺瘤、
# 把细胞系的 "control" 判成正常组织（两处均为实测踩坑）。
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from cohort_hunt2 import classify as _classify  # noqa
except Exception:  # pragma: no cover
    _classify = None


def classify3(title, source="", chars=""):
    """三字段分类；有 cohort_hunt2 的验证版就用它，否则退回朴素法。"""
    if _classify is not None:
        return _classify(title, source, chars)[0]
    s = (title + " " + source + " " + chars).lower()
    if any(k in s for k in NORM):
        return "normal"
    if ("carcinoma" in s or "adenocarcinoma" in s or "cancer" in s or "tumor" in s
            or "tumour" in s or re.search(r"\bcrc\b", s) or "malignan" in s):
        return "carcinoma"
    if "adenoma" in s or "adenomatous" in s or "polyp" in s:
        return "adenoma"
    return "other"



def parse_series(gse):
    raw = fetch(matrix_url(gse))
    if not raw:
        return None
    try:
        txt = gzip.decompress(raw).decode("utf-8", "replace")
    except Exception:
        return None
    lines = txt.split("\n")
    titles, srcs, gpl = [], [], None
    chars = []
    tbl_at = None
    for i, ln in enumerate(lines):
        if ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_source_name_ch1"):
            srcs = [x.strip().strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip().strip('"') for x in ln.split("\t")[1:]])
        elif ln.startswith("!Sample_platform_id"):
            gpl = ln.split("\t")[1].strip().strip('"')
        elif ln.startswith("!series_matrix_table_begin"):
            tbl_at = i + 1
            break
    if tbl_at is None:
        return None
    ns = max(len(titles), len(srcs))
    char_by_s = [""] * ns
    for row in chars:
        for j in range(min(ns, len(row))):
            char_by_s[j] += " " + row[j]
    labels = [classify3(titles[i] if i < len(titles) else "",
                        srcs[i] if i < len(srcs) else "",
                        char_by_s[i])
              for i in range(ns)]
    rdr = csv.reader(lines[tbl_at:], delimiter="\t")
    hdr = next(rdr, None)
    if not hdr:
        return None
    expr = {}
    for row in rdr:
        if not row or row[0].startswith("!"):
            break
        pid = row[0].strip().strip('"')
        vals = []
        for v in row[1:]:
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(None)
        if len(vals) >= 6:
            expr[pid] = vals
    return {"labels": labels, "expr": expr, "gpl": gpl, "titles": titles}


def load_entrez_map(row_ids):
    """
    备用注释路径：当平台既无 .annot.gz 也无预生成映射表时，看行 ID 是不是 Entrez 形式。
    实测：GSE117606 的平台 GPL25373 行 ID 形如 `100009676_at` / `10000_at`
    （已按 Entrez Gene ID 重新注释），去掉 `_at` 即 Entrez ID。
    映射来源：geo_cache/Homo_sapiens.gene_info.gz（NCBI，5.2 MB，制表符分隔，
              列：tax_id / GeneID / Symbol / ...），比 55.8 MB 的平台 SOFT 划算得多。
    """
    ids = [str(x) for x in row_ids]
    hit = [x for x in ids[:2000] if re.match(r"^\d+_at$", x)]
    if len(hit) < 50:
        return {}
    gi = os.path.join(CACHE, "Homo_sapiens.gene_info.gz")
    if not os.path.exists(gi):
        print("  [annot] 需要 Entrez 映射但缺少 Homo_sapiens.gene_info.gz", flush=True)
        return {}
    out = {}
    with gzip.open(gi, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) > 2 and p[0] == "9606":     # 人类
                out[p[1]] = p[2]
    mapped = {}
    for x in ids:
        m = re.match(r"^(\d+)_at$", x)
        if m and m.group(1) in out:
            mapped[x] = out[m.group(1)]
    print("  [annot] Entrez 映射：%d / %d 行 ID 已解析为基因符号" % (len(mapped), len(ids)), flush=True)
    return mapped


def load_annot(gpl):
    """
    探针 → 基因符号。两级来源：
      ① 预生成的 geo_cache/GPLxxx_map.tsv（由 gpl_soft_map.py 从 family.soft.gz 抽出）——
         用于没有官方 .annot.gz 的平台（如 Agilent GPL25373）。
      ② 官方 GPLxxx.annot.gz（只有部分平台有；实测 GPL570 / GPL4133 有）。
    列名兼容：'Gene symbol' 与 'GENE_SYMBOL' 两种写法都存在（Agilent 是后者）。
    """
    if not gpl:
        return {}
    tsv = os.path.join(CACHE, gpl + "_map.tsv")
    if os.path.exists(tsv) and os.path.getsize(tsv) > 1000:
        out = {}
        for ln in open(tsv, encoding="utf-8"):
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0]:
                out[p[0]] = p[1]
        if out:
            print("  [annot] 使用预生成映射表 %s（%d 条）" % (os.path.basename(tsv), len(out)), flush=True)
            return out
    raw = fetch(gpl_url(gpl))
    if not raw:
        return {}
    try:
        txt = gzip.decompress(raw).decode("utf-8", "replace")
    except Exception:
        return {}
    out = {}
    started = False
    hidx, iidx = 1, 0
    for ln in txt.split("\n"):
        if ln.startswith("!platform_table_begin"):
            started = True
            continue
        if ln.startswith("!platform_table_end"):
            break
        if not started:
            continue
        parts = ln.split("\t")
        if not out and any(p.strip().lower().replace("_", " ") == "gene symbol" for p in parts):
            for j, p in enumerate(parts):
                if p.strip().lower().replace("_", " ") == "gene symbol":
                    hidx = j
                if p.strip().lower() == "id":
                    iidx = j
            continue
        if len(parts) <= max(hidx, iidx):
            continue
        sym = parts[hidx].strip().strip('"')
        if not sym or sym in ("", "---"):
            continue
        out[parts[iidx].strip().strip('"')] = sym.split("///")[0].strip()
    return out


def auc(pos, neg):
    if not pos or not neg:
        return None
    allv = sorted([(v, 0) for v in neg] + [(v, 1) for v in pos])
    # rank with ties
    ranks = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1
        ranks[allv[i][0]] = r
        i = j + 1
    n1, n2 = len(pos), len(neg)
    R1 = sum(ranks[v] for v in pos)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    return U1 / (n1 * n2)


def fmt_auc(a, pos="癌", neg="腺瘤"):
    """
    ★ 显示修正：AUC<0.5 不是「没有区分力」，而是「反向完美区分」。
      实测踩坑：下调分量在 LCM 队列得 AUC=0.000（完美反向），旧代码 `if aa else "—"`
      把 0.0 当成假值印成了「—」，差点把最强的一条证据读丢。统一按 max(a,1-a) 报「区分力」，
      并显式标出方向。
    """
    if a is None:
        return "—"
    d = max(a, 1 - a)
    if abs(a - 0.5) < 1e-9:
        return "%.3f (无方向差异)" % d
    return "%.3f (%s > %s)" % (d, pos if a > 0.5 else neg, neg if a > 0.5 else pos)


SIG = list(csv.DictReader(open(SIGFILE, encoding="utf-8-sig")))
UP = [r["gene"] for r in SIG if float(r["mean_logFC"]) > 0]
DN = [r["gene"] for r in SIG if float(r["mean_logFC"]) < 0]
print("签名：%d 基因（上调 %d / 下调 %d）" % (len(SIG), len(UP), len(DN)))

PANELS = {
    "间质·激活/CAF": "FAP THBS2 COL10A1 COL11A1 INHBA COMP TNC LRRC15 CTHRC1 SULF1 SULF2 GREM1 WNT5A POSTN ASPN MMP11 MMP13 MMP1 LOXL2 CDH11 PXDN ADAM12 SEMA3C NID2 ADAMTS2 ISLR EMILIN1 COL8A1 SPARC THY1 SERPINE1 TIMP1 MMP3 MMP14 LOX PLOD2 BGN FNDC1 COL16A1 THBS1 CTGF ANGPTL4 VEGFA",
    "间质·基线": "DCN LUM PDGFRA FBLN1 FBLN2 OGN PRELP FMOD C7 CFH SERPINF1 PCOLCE MFAP4 MFAP5 ECM1 PTGDS SOD3 MMP2 IGFBP7 CCN1 CCN2 IGFBP5 GAS6 LAMA4 LAMB1 ITGA1 CAV1 CAV2 RARRES2 NNMT NBL1",
    "上皮·分化/功能": "MUC2 MUC4 MUC12 MUC13 MUC20 TFF1 TFF2 TFF3 SI OLFM4 PIGR ZG16 CA1 CA2 CEACAM7 GUCA2A GUCA2B AQP8 SLC26A3 SLC51B HPGD CD177 CLCA4 PLAC8 FABP1 KRT20 MEP1A ANPEP APOA4 CDX2 VIL1 REG4 SPINK4 FCGBP AGR2 ITLN1 HNF4A KRT23 ADH1C SLC12A2 MUC1",
    "上皮·基线/结构": "EPCAM KRT8 KRT18 KRT19 CDH1 CLDN3 CLDN4 CLDN7 DSP JUP PKP3 ITGB4 ITGA6 SFN RAB25 PRSS8 GRHL2 ELF3 EHF ESRP1 ESRP2 DSC2 PKP2 OCLN",
    "髓系·激活": "SPP1 CCL18 CCL2 CCL3 CCL4 CCL20 APOE TREM2 MMP9 MMP12 SIGLEC1 CXCL8 IL10 IDO1 CD274 PDCD1LG2 LGALS3 S100A8 S100A9 S100A12 FCN1 C5AR1 MNDA FCGR3B MMP8 TCN1 CCL19",
}
PANELS = {k: set(v.split()) for k, v in PANELS.items()}

results = {}
# ★ 不带参数时必须**默认跑这四个队列**，不能空转。
#   旧写法 `for gse in sys.argv[1:]` 在无参数时一次都不循环，末尾却照样 dump ——
#   于是 `python sig_validate.py` 把 cohort_hunt/sig_validation_results.json 覆盖成 `{}`，
#   而这个文件正是 gse117606_paired.py 要备份、manuscript_facts.py → build_suppl_tables.py
#   要读的四队列 AUC 来源（README 的第一个 smoke test 就会踩到）。清空后不报错，
#   只是四个队列的 AUC 从此消失。
COHORTS_4 = ["GSE4183", "GSE117606", "GSE128435", "GSE15960"]
# ★ 主循环只在**直接运行**时执行。three_layer_check / gse117606_paired / lcm_diag 都是用
#   exec_module 载入本文件「借用函数」，导入即跑会让它们各自重复下载四个队列。
_TODO = (sys.argv[1:] or COHORTS_4) if __name__ == "__main__" else []
for gse in _TODO:
    print("\n" + "=" * 108)
    print("验证队列 %s" % gse)
    print("=" * 108)
    d = parse_series(gse)
    if d is None:
        print("  !! 无法解析 series matrix（可能非芯片平台）")
        results[gse] = {"ok": False}
        continue
    cnt = Counter(d["labels"])
    print("  平台 %s | 样本 %d | 腺瘤 %d / 癌 %d / 正常 %d / 其他 %d"
          % (d["gpl"], len(d["labels"]), cnt.get("adenoma", 0), cnt.get("carcinoma", 0),
             cnt.get("normal", 0), cnt.get("other", 0)))
    # 阈值放到 4：GSE15960 只有 5 腺瘤/6 癌，但它是**激光显微切割(LCM)上皮细胞**队列，
    # 是验证「上调分量来自间质、下调分量是上皮内在」的天然对照，不能因为例数少就丢掉。
    MINN = int(os.environ.get("SIGVAL_MIN", "4"))
    if cnt.get("adenoma", 0) < MINN or cnt.get("carcinoma", 0) < MINN:
        print("  !! 腺瘤或癌例数不足(%d)，跳过" % MINN)
        results[gse] = {"ok": False, "counts": dict(cnt)}
        continue
    ann = load_annot(d["gpl"])
    if not ann:
        # 备用路径：行 ID 为 Entrez 形式（如 GSE117606 的 `10000_at`）时用 gene_info 映射
        ann = load_entrez_map(d["expr"].keys())
    print("  注释探针 %d 条" % len(ann))
    if not ann:
        print("  !! 无法取得平台注释")
        results[gse] = {"ok": False, "counts": dict(cnt)}
        continue

    # probe → gene（多探针取跨样本均值最大者）
    gene_rows = {}
    for pid, vals in d["expr"].items():
        g = ann.get(pid)
        if not g:
            continue
        good = [v for v in vals if v is not None]
        if len(good) < len(vals) * 0.8:
            continue
        m = statistics.mean(good)
        if g not in gene_rows or m > gene_rows[g][0]:
            gene_rows[g] = (m, vals)
    print("  可比对基因 %d 个" % len(gene_rows))

    # log2 变换（若非 log 尺度）
    sample_med = []
    for i in range(len(d["labels"])):
        v = [gene_rows[g][1][i] for g in gene_rows
             if i < len(gene_rows[g][1]) and gene_rows[g][1][i] is not None]
        if v:
            sample_med.append(statistics.median(v))
    if sample_med and statistics.median(sample_med) > 30:
        print("  检测到线性尺度 → 做 log2(x+1)")
        for g in gene_rows:
            gene_rows[g] = (gene_rows[g][0], [math.log2(v + 1) if v is not None else None
                                              for v in gene_rows[g][1]])

    # 基因级 z（队列内）
    ng = len(d["labels"])


    def zrow(vals):
        good = [v for v in vals[:ng] if v is not None]
        if len(good) < 3:
            return None
        m = statistics.mean(good)
        sd = statistics.pstdev(good) or 1.0
        return [(v - m) / sd if v is not None else 0.0 for v in vals[:ng]]


    zmap = {}
    for g in gene_rows:
        z = zrow(gene_rows[g][1])
        if z:
            zmap[g] = z

    idx_ad = [i for i, l in enumerate(d["labels"]) if l == "adenoma"][:ng]
    idx_ca = [i for i, l in enumerate(d["labels"]) if l == "carcinoma"][:ng]

    # 签名分数
    def score(i):
        up = [zmap[g][i] for g in UP if g in zmap]
        dn = [zmap[g][i] for g in DN if g in zmap]
        if not up or not dn:
            return None
        return statistics.mean(up) - statistics.mean(dn)

    s_ad = [score(i) for i in idx_ad]
    s_ca = [score(i) for i in idx_ca]
    s_ad = [x for x in s_ad if x is not None]
    s_ca = [x for x in s_ca if x is not None]
    a = auc(s_ca, s_ad)
    print("\n  ★ 签名分（上调z均值 − 下调z均值）")
    print("     腺瘤 n=%d 中位 %+.3f   癌 n=%d 中位 %+.3f"
          % (len(s_ad), statistics.median(s_ad) if s_ad else float('nan'),
             len(s_ca), statistics.median(s_ca) if s_ca else float('nan')))
    print("     癌 vs 腺瘤 AUC = %s" % fmt_auc(a))

    # ---- 分量 AUC：把签名拆成「上调分量(间质/髓系激活)」与「下调分量(上皮功能关闭)」----
    # 这一拆分才是可发表的点：若上调分量单独就有高 AUC，说明「间质激活」足以区分腺瘤与癌，
    # 这正是 mIF 面板该染什么的数据依据（而全签名 AUC 高是必然的，不构成发现）。
    def score_part(i, genes):
        zz = [zmap[g][i] for g in genes if g in zmap]
        return statistics.mean(zz) if zz else None

    comp = {}
    for tag, genes in (("上调分量(间质/髓系激活)", UP), ("下调分量(上皮功能关闭)", DN)):
        p_ad = [score_part(i, genes) for i in idx_ad]
        p_ca = [score_part(i, genes) for i in idx_ca]
        p_ad = [x for x in p_ad if x is not None]
        p_ca = [x for x in p_ca if x is not None]
        aa = auc(p_ca, p_ad)
        comp[tag] = aa
        print("     %-22s n可用基因=%3d  AUC = %-22s (腺瘤中位 %+.3f / 癌中位 %+.3f)"
              % (tag, len([g for g in genes if g in zmap]), fmt_auc(aa),
                 statistics.median(p_ad) if p_ad else float('nan'),
                 statistics.median(p_ca) if p_ca else float('nan')))

    # ---- 面板级 AUC：每个细胞面板单独能否区分腺瘤与癌 ----
    print("\n  ★ 各细胞面板单独 AUC（0.5=无区分力）")
    panel_auc = {}
    for k, genes in PANELS.items():
        p_ad = [score_part(i, genes) for i in idx_ad]
        p_ca = [score_part(i, genes) for i in idx_ca]
        p_ad = [x for x in p_ad if x is not None]
        p_ca = [x for x in p_ca if x is not None]
        if len(p_ad) < 3 or len(p_ca) < 3:
            continue
        aa = auc(p_ca, p_ad)
        panel_auc[k] = aa
        print("     %-14s n=%3d  AUC = %s" % (k, len([g for g in genes if g in zmap]), fmt_auc(aa)))

    # ---- ★ 关键对照：腺瘤 vs 正常 ----
    # 目的：判断签名到底在「正常→腺瘤」就已经变了，还是「腺瘤→癌」才变。
    # 若上调分量在正常→腺瘤就已升高（AUC 高），说明间质激活是**腺瘤期事件**，
    # 那么腺瘤→癌的进一步升高只可能是「间质占比继续增加」（组成效应），
    # 这与 Front Oncol 2022 的蛋白组学结论（腺瘤与癌间质可比）互相印证。
    idx_no = [i for i, l in enumerate(d["labels"]) if l == "normal"][:ng]
    an_auc = {}
    if len(idx_no) >= 4:
        print("\n  ★ 对照：腺瘤 vs 正常（判断签名在序贯的哪一步开始变）")
        s_no = [score(i) for i in idx_no]
        s_no = [x for x in s_no if x is not None]
        a2 = auc(s_ad, s_no)
        an_auc["signature"] = a2
        print("     全签名           A 腺瘤 vs 正常 AUC = %s" % fmt_auc(a2, "腺瘤", "正常"))
        for tag, genes in (("上调分量(间质/髓系激活)", UP), ("下调分量(上皮功能关闭)", DN)):
            p_no = [score_part(i, genes) for i in idx_no]
            p_no = [x for x in p_no if x is not None]
            # 注意：必须在这里重新算腺瘤分量分，不能沿用上面面板循环残留的 p_ad
            p_ad_ = [score_part(i, genes) for i in idx_ad]
            p_ad_ = [x for x in p_ad_ if x is not None]
            aa = auc(p_ad_, p_no)
            an_auc[tag] = aa
            print("     %-22s AUC = %-22s (正常中位 %+.3f / 腺瘤中位 %+.3f)"
                  % (tag, fmt_auc(aa, "腺瘤", "正常"),
                     statistics.median(p_no) if p_no else float('nan'),
                     statistics.median(p_ad_) if p_ad_ else float('nan')))

    # 逐基因方向符合
    hit = tot = 0
    panel_fc = {k: [] for k in PANELS}
    for r in SIG:
        g = r["gene"]
        if g not in zmap:
            continue
        tot += 1
        exp = 1 if float(r["mean_logFC"]) > 0 else -1
        obs = 1 if statistics.mean([zmap[g][i] for i in idx_ca]) > \
                   statistics.mean([zmap[g][i] for i in idx_ad]) else -1
        if exp == obs:
            hit += 1
    for k, s in PANELS.items():
        for g in s:
            if g in zmap:
                fc = (statistics.mean([zmap[g][i] for i in idx_ca]) -
                      statistics.mean([zmap[g][i] for i in idx_ad]))
                panel_fc[k].append((g, fc))
    print("\n  逐基因方向符合率：%d / %d = %.1f%%" % (hit, tot, 100.0 * hit / max(tot, 1)))
    print("\n  各面板中位 Δ（癌 − 腺瘤，单位=队列内 z 差；符号应复现构建期结论）")
    for k in ["间质·激活/CAF", "髓系·激活", "间质·基线", "上皮·分化/功能", "上皮·基线/结构"]:
        L = panel_fc[k]
        if not L:
            print("     %-14s 0 命中" % k)
            continue
        med = statistics.median([x[1] for x in L])
        ups = sum(1 for x in L if x[1] > 0)
        print("     %-14s %3d 命中 | 上调 %2d | 中位 Δ %+6.3f | %s"
              % (k, len(L), ups, med,
                 ", ".join("%s(%+.2f)" % x for x in sorted(L, key=lambda y: -abs(y[1]))[:5])))
    results[gse] = {"ok": True, "gpl": d["gpl"], "counts": dict(cnt),
                    "auc": a, "component_auc": comp, "panel_auc": panel_auc,
                    "adenoma_vs_normal_auc": an_auc,
                    "gene_dir_concordance": hit / max(tot, 1), "n_tested": tot,
                    "score_adenoma_median": statistics.median(s_ad) if s_ad else None,
                    "score_carcinoma_median": statistics.median(s_ca) if s_ca else None,
                    "panel_median_delta": {k: (statistics.median([x[1] for x in v]) if v else None)
                                           for k, v in panel_fc.items()}}

OUT = os.path.join(BASE, "cohort_hunt", "sig_validation_results.json")
if __name__ != "__main__":
    pass                      # 被 import：不写盘，留给调用方自己的产物
elif not results:
    # 一个队列都没算成（缺数据/注释）：**绝不覆盖**已有结果，否则下游静默丢掉四队列 AUC
    print("\n[skip] 一个队列都没算成，保留既有 %s 不动" % os.path.basename(OUT))
else:
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(results, open(OUT, "w", encoding="utf-8"),
              ensure_ascii=False, indent=2, default=str)
    print("\n[saved] %s（%d 个队列）" % (os.path.basename(OUT), len(results)))
print("DONE")
