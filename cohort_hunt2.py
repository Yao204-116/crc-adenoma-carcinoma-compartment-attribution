# -*- coding: utf-8 -*-
"""
搜寻「第 5 个独立队列」——即同时含腺瘤(adenoma/polyp) 与 癌(carcinoma) 且不在训练集里的 GSE。

设计要点（吸取上一版 cohort_hunt.py 的教训）：
  * 上一版用 acc.cgi?targ=all 拿样本元数据 -> 巨慢 + IncompleteRead，跑 31 分钟被终止。
  * 实测 acc.cgi?targ=gsm&form=text&view=brief 一次请求返回全部 GSM 记录（3~6 秒），
    含 !Sample_title / !Sample_source_name_ch1 / !Sample_characteristics_ch1 / !Sample_platform_id
    -> 这是唯一需要的接口。
  * 关键原则（本项目反复踩坑）：**GEO 命中 ≠ 内容符合**。
    绝不能凭 Series 标题判断有没有腺瘤样本，必须逐条 GSM 分类计数。
  * 分类用带词边界的正则：`\\badenoma\\b` 不会误匹配 adenocarcinoma（因为 'c' 是词字符）。

输出：
  cohort_hunt2/gsm_cache/GSE*.txt  原始缓存
  cohort_hunt2/cohort_scan.json    全部候选的逐样本分类结果
  cohort_hunt2/cohort_usable.json  只保留 腺瘤>=5 且 癌>=5 的可用队列（按腺瘤数排序）
"""
import concurrent.futures as cf
import json, os, re, sys, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cohort_hunt2")
CACHE = os.path.join(OUT, "gsm_cache")
os.makedirs(CACHE, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0"}
EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# 已经用过的队列（训练集），必须排除，否则验证就是循环论证
TRAIN = {"GSE41258", "GSE37364", "GSE20916", "GSE41657"}

# ---- 标签正则 ----------------------------------------------------------
# 关键：先把文本里的非字母数字统一成空格，这样 "colon_normal_1081" -> "colon normal 1081"、
# "CRC_CD1" -> "crc cd1"，词边界才能真正生效（下划线在正则里是词字符，会破坏 \b）。
def norm(t):
    return re.sub(r"[^a-z0-9]+", " ", (t or "").lower())


RE_CELL = re.compile(r"cell ?line|caco|hct ?116|sw ?480|sw ?620|ht ?29|\brko\b|\bdld ?1|"
                     r"lovo|colo ?205|colo ?320|hke3|lim ?1215|fet\b|organoid", re.I)
# 明确表示「非肿瘤 / 瘤旁」的措辞，优先级最高（必须挡在 carcinoma 前面）
RE_NONNEO = re.compile(r"adjacent|non ?tumou?r|nontumou?r|non ?neoplastic|uninvolved|"
                       r"peri ?tumou?r|normal|healthy|distant|benign", re.I)
RE_ADEN = re.compile(r"adenoma|adenomatous|polyp", re.I)   # 不含 adenocarcinoma（无 'adenoma' 子串）
RE_CARC = re.compile(r"adenocarcinoma|carcinoma|\bcrc\b|malignan|colon cancer|colorectal cancer|"
                     r"rectal cancer|tumou?r|cancer|neoplas", re.I)

# 注意：**取消了 \bcontrol\b** —— 它把细胞系的 "cell line, control" 误判成正常组织（实测 GSE32323 多出 4 例）


def label_of(field):
    """对单个字段判定；无信号返回 None。字段内优先级：细胞系 > 非肿瘤 > 癌 > 腺瘤。"""
    s = norm(field)
    if not s.strip():
        return None
    if RE_CELL.search(s):
        return "cellline"
    if RE_NONNEO.search(s):
        return "normal"
    if RE_CARC.search(s):
        return "carcinoma"
    if RE_ADEN.search(s):
        return "adenoma"
    return None


def fetch(url, timeout=90):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) \
        .read().decode("utf-8", "replace")


def get_json(url, timeout=60, tries=3):
    for i in range(tries):
        try:
            return json.loads(fetch(url, timeout))
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(1.5 * (i + 1))


# ---------------- Phase 1: 发现候选 GSE ----------------
QUERIES = [
    "(colorectal OR colon OR rectal) AND (adenoma OR adenomatous OR polyp) AND "
    "(carcinoma OR adenocarcinoma OR cancer) AND Homo sapiens[Organism] AND gse[Entry Type]",
    "adenoma AND carcinoma AND (colon OR colorectal OR rectal) AND gse[Entry Type]",
    '"adenoma-carcinoma sequence" AND gse[Entry Type]',
    "(adenoma OR polyp) AND (carcinoma OR cancer) AND colon AND expression profiling[DataSet Type] AND gse[Entry Type]",
    "(adenoma OR adenomatous polyp) AND (adenocarcinoma) AND Homo sapiens[Organism] AND gse[Entry Type]",
    "colorectal neoplasm AND (normal OR adenoma) AND (carcinoma OR adenocarcinoma) AND gse[Entry Type]",
]

# 额外手工种子（已知可能含腺瘤+癌，用于兜底，仍需逐样本核实）
SEED = [
    "GSE4183", "GSE8671", "GSE18105", "GSE10714", "GSE23878", "GSE9348",
    "GSE4107", "GSE22598", "GSE50118", "GSE24514", "GSE77953", "GSE100179",
    "GSE113513", "GSE117606", "GSE136430", "GSE44076", "GSE50118", "GSE37250",
    "GSE26906", "GSE15960", "GSE30292", "GSE142279", "GSE156720",
]


def discover():
    """
    发现候选 GSE。
    ★ 实测教训：db=gds 的 esummary 耗时对批大小**非线性**爆炸 ——
      200 个 id = 98.8 秒（0.49 s/id），而 40 个 id 只要 4.6 秒（0.12 s/id）。
      所以必须用小批（40），大批会让整个发现阶段卡死（上一版就是这么超时的）。
    """
    found = {}

    def one_query(q):
        out = {}
        try:
            u = "%s/esearch.fcgi?db=gds&retmode=json&retmax=500&term=%s" % (EUT, urllib.parse.quote(q))
            r = get_json(u, 60)
            ids = r["esearchresult"]["idlist"]
            if not ids:
                return out
            for i in range(0, len(ids), 40):
                batch = ids[i:i + 40]
                s = get_json("%s/esummary.fcgi?db=gds&retmode=json&id=%s"
                             % (EUT, ",".join(batch)), 90)["result"]
                for k in s["uids"]:
                    d = s[k]
                    acc = d.get("accession", "")
                    if acc.startswith("GSE"):
                        out.setdefault(acc, {
                            "title": d.get("title", ""),
                            "n_samples": int(d.get("n_samples") or 0),
                            "gpl": d.get("gpl", ""),
                            "pdat": d.get("pdat", ""),
                        })
                time.sleep(0.4)   # 尊重 E-utilities 速率限制
        except Exception as e:
            print("  [warn] query failed: %s -> %s" % (q[:46], type(e).__name__), flush=True)
        return out

    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        for out in ex.map(one_query, QUERIES):
            for k, v in out.items():
                found.setdefault(k, v)
            print("  candidates=%d" % len(found), flush=True)

    for acc in SEED:
        found.setdefault(acc, {"title": "(seed)", "n_samples": 0, "gpl": "", "pdat": ""})
    return found


# ---------------- Phase 2: 逐样本分类 ----------------
def classify(title, source, chars):
    """
    优先级 + 回退：title 最能反映该样本的真实身份；title 无信号才看 source，最后看 characteristics。
    理由（实测教训）：很多 GSE 的 !Sample_characteristics_ch1 里含**研究级字段**
    （如 GSE8671 所有样本都带腺瘤相关的统一标注），直接拼接会把正常黏膜误判成腺瘤。
    返回 (label, 命中字段)。
    """
    for fld, name in ((title, "title"), (source, "source"), (chars, "chars")):
        lab = label_of(fld)
        if lab:
            return lab, name
    return "other", ""


def scan_one(acc):
    """拉 GSE 的全部 GSM 记录并逐样本分类。返回 dict。"""
    cf_ = os.path.join(CACHE, acc + ".txt")
    try:
        if os.path.exists(cf_) and os.path.getsize(cf_) > 200:
            t = open(cf_, encoding="utf-8", errors="replace").read()
        else:
            t = fetch("%s?acc=%s&targ=gsm&form=text&view=brief" % (ACC, acc), timeout=120)
            open(cf_, "w", encoding="utf-8").write(t)
    except Exception as e:
        return {"acc": acc, "error": "%s: %s" % (type(e).__name__, str(e)[:100])}

    blocks = re.split(r"\n\^SAMPLE\s*=\s*", t)
    m = re.findall(r"!Sample_platform_id\s*=\s*(\S+)", t)
    plat = m[0] if m else ""

    labels, samples = [], []
    for b in blocks[1:]:
        gsm = b.split("\n", 1)[0].strip()
        title = source = ""
        chars = []
        for ln in b.split("\n"):
            if ln.startswith("!Sample_title"):
                title = ln.split("=", 1)[1]
            elif ln.startswith("!Sample_source_name"):
                source = ln.split("=", 1)[1]
            elif ln.startswith("!Sample_characteristics"):
                chars.append(ln.split("=", 1)[1])
        chars = " ".join(chars)
        lab, via = classify(title.replace('"', " "), source.replace('"', " "), chars.replace('"', " "))
        labels.append(lab)
        samples.append({"gsm": gsm, "label": lab, "via": via,
                        "title": title.strip()[:90],
                        "chars": chars.strip()[:110]})

    cnt = {k: labels.count(k) for k in ("normal", "adenoma", "carcinoma", "cellline", "other")}
    return {"acc": acc, "platform": plat, "n": len(labels), "counts": cnt, "samples": samples}


def main():
    print("=== Phase 1: 发现候选 ===", flush=True)
    cand = discover()
    cands = sorted(cand.keys())
    print("候选总数: %d" % len(cands), flush=True)

    todo = [a for a in cands if a not in TRAIN]
    print("排除训练集后待扫描: %d" % len(todo), flush=True)

    results = []
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(scan_one, a): a for a in todo}
        done = 0
        for f in cf.as_completed(futs):
            done += 1
            r = f.result()
            results.append(r)
            if r.get("error"):
                tag = "ERR " + r["error"][:60]
            else:
                c = r["counts"]
                tag = "AD=%d CA=%d NO=%d CL=%d" % (c["adenoma"], c["carcinoma"], c["normal"], c["cellline"])
            print("[%3d/%d] %-10s %-34s (%.0fs)" % (done, len(todo), r["acc"], tag, time.time() - t0), flush=True)

    for r in results:
        if not r.get("error"):
            r["title"] = cand.get(r["acc"], {}).get("title", "")
            r["pdat"] = cand.get(r["acc"], {}).get("pdat", "")

    ok = [r for r in results if not r.get("error")]
    usable = [r for r in ok if r["counts"]["adenoma"] >= 5 and r["counts"]["carcinoma"] >= 5]
    usable.sort(key=lambda r: -r["counts"]["adenoma"])

    json.dump({"candidates": cand, "results": results},
              open(os.path.join(OUT, "cohort_scan.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    json.dump(usable, open(os.path.join(OUT, "cohort_usable.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)

    print("\n=== 扫描完成 用时 %.0f 分钟 ===" % ((time.time() - t0) / 60), flush=True)
    print("成功扫描 %d 个 GSE；其中 腺瘤>=5 且 癌>=5 的可用队列 %d 个：" % (len(ok), len(usable)), flush=True)
    print("%-11s %-6s %-5s %-5s %-5s %-5s %-5s  %s" % ("GSE", "平台", "腺瘤", "癌", "正常", "细胞系", "其他", "标题"), flush=True)
    for r in usable:
        c = r["counts"]
        print("%-11s %-6s %-5d %-5d %-5d %-5d %-5d  %s" % (
            r["acc"], r.get("platform", "?"), c["adenoma"], c["carcinoma"], c["normal"],
            c["cellline"], c["other"], (r.get("title") or "")[:58]), flush=True)

    # 额外：只有癌+正常、或只有腺瘤+正常的大队列也值得记一笔
    print("\n--- 参考：样本量>=40 的其他结肠队列（无腺瘤或未分类成功） ---", flush=True)
    for r in sorted(ok, key=lambda x: -x["n"]):
        c = r["counts"]
        if r["n"] >= 40 and r not in usable:
            print("%-11s n=%-4d AD=%-3d CA=%-3d NO=%-3d CL=%-3d OT=%-3d %s" % (
                r["acc"], r["n"], c["adenoma"], c["carcinoma"], c["normal"], c["cellline"], c["other"],
                (r.get("title") or "")[:46]), flush=True)


if __name__ == "__main__":
    main()
