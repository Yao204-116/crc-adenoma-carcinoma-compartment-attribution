# -*- coding: utf-8 -*-
"""GSE303070 全量元数据核查：是否与「腺瘤期分化程序内部结构」正面撞车。

上轮线索：n=352 空间图谱，摘要提到 "stromal niches on tumour differentiation"。
本脚本拉：series 摘要全文、GSM 构成、平台、关联 PMID。
"""
import re, gzip, io, json, time, urllib.parse, urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research; mailto:test@example.com)"}
EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"


def _get(u, t=90, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=t) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 + 3 * i)


out = {}

# 1) series 元数据（brief）
u = BASE + urllib.parse.urlencode({"acc": "GSE303070", "targ": "self", "form": "text", "view": "brief"})
txt = _get(u).decode("utf-8", "replace")
out["series_brief"] = txt
print("=== SERIES BRIEF ===")
print(txt[:3000])

# 2) 完整摘要
u2 = BASE + urllib.parse.urlencode({"acc": "GSE303070", "targ": "self", "form": "text", "view": "quick"})
txt2 = _get(u2).decode("utf-8", "replace")
m = re.search(r"!Series_summary = (.*?)(?=\n!Series_)", txt2, re.S)
out["summary"] = m.group(1) if m else ""
print("\n=== SUMMARY ===")
print(out["summary"][:4000])

m = re.search(r"!Series_overall_design = (.*?)(?=\n!Series_)", txt2, re.S)
out["design"] = m.group(1) if m else ""
print("\n=== OVERALL DESIGN ===")
print(out["design"][:3000])

for k in ("Series_title", "Series_platform_id", "Series_pubmed_id", "Series_sample_id",
          "Series_type", "Series_relation"):
    for mm in re.finditer(r"!%s = (.*)" % k, txt2):
        print("%s = %s" % (k, mm.group(1)))

# 3) GSM 构成
u3 = BASE + urllib.parse.urlencode({"acc": "GSE303070", "targ": "gsm", "form": "text", "view": "brief"})
txt3 = _get(u3).decode("utf-8", "replace")
titles = re.findall(r"!Sample_title = (.*)", txt3)
chars = re.findall(r"!Sample_characteristics_ch1 = (.*)", txt3)
out["n_gsm"] = len(titles)
out["titles"] = titles[:60]
out["chars_sample"] = chars[:80]
print("\n=== GSM n=%d ===" % len(titles))
for t in titles[:40]:
    print("  ", t)
from collections import Counter
print("\n-- chars 词频 --")
for k, v in Counter(chars).most_common(40):
    print("  %4d  %s" % (v, k))

# 4) pubmed
pmids = re.findall(r"!Series_pubmed_id = (\d+)", txt2)
out["pmids"] = pmids
if pmids:
    u4 = EUT + "efetch.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "text"})
    abs_txt = _get(u4).decode("utf-8", "replace")
    out["pubmed"] = abs_txt
    print("\n=== PUBMED ===")
    print(abs_txt[:5000])

# 5) supplementary 体积
sups = re.findall(r"!Series_supplementary_file = (.*)", txt2)
out["series_suppl"] = sups
print("\n=== SERIES SUPPL ===")
for s in sups:
    print("  ", s)

json.dump(out, open("gse303070_probe.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("\nsaved gse303070_probe.json")
