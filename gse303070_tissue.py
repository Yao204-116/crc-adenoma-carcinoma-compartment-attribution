# -*- coding: utf-8 -*-
"""核实 GSE303070 的样本组织构成（tumour vs normal，有无腺瘤）。"""
import urllib.parse, urllib.request, time, re, json, sys

UA = {"User-Agent": "Mozilla/5.0 (research)"}


def get(u, t=120, tries=5):
    for i in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=t).read().decode("utf-8", "replace")
        except Exception as e:
            print("retry", i, repr(e)[:90], flush=True)
            time.sleep(4)
    return ""


# 1) 全样本 brief 文本
u = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?" + urllib.parse.urlencode(
    {"acc": "GSE303070", "targ": "gsm", "form": "text", "view": "brief"})
txt = get(u)
open("gse303070_gsm_brief.txt", "w", encoding="utf-8").write(txt)
print("brief bytes:", len(txt), flush=True)

# 2) 逐 GSM 解析 Sample_title / tissue
blocks = re.split(r"\^SAMPLE = ", txt)[1:]
print("n_GSM parsed:", len(blocks), flush=True)

rows = []
for b in blocks:
    gsm = b.split("\n", 1)[0].strip()
    ti = re.search(r"!Sample_title = (.+)", b)
    ch = re.findall(r"!Sample_characteristics_ch1 = (.+)", b)
    tissue = [c for c in ch if c.lower().startswith("tissue:")]
    donor = [c for c in ch if c.lower().startswith("donor:")]
    cells = [c for c in ch if "total cells" in c.lower()]
    rows.append({
        "gsm": gsm,
        "title": ti.group(1).strip() if ti else "",
        "tissue": tissue[0].split(":", 1)[1].strip() if tissue else "",
        "donor": donor[0].split(":", 1)[1].strip() if donor else "",
        "cells": int(re.search(r"(\d+)", cells[0]).group(1)) if cells else 0,
    })

# 3) 汇总：组织类型 / 病例(field) / FOV
from collections import Counter, defaultdict
tc = Counter(r["tissue"] for r in rows)
print("\n=== tissue 字段分布 ===", flush=True)
for k, v in tc.most_common():
    print("  %-40s %d" % (k, v), flush=True)

# title 前缀（FOV 前的字段）→ 病例单位
def case_key(t):
    m = re.match(r"([^-]+)", t)
    return m.group(1) if m else t
cc = Counter(case_key(r["title"]) for r in rows)
print("\n=== title 首段（病例/切片标识）%d 个唯一值 ===" % len(cc), flush=True)
for k, v in sorted(cc.items()):
    print("  %-30s %d FOV" % (k, v), flush=True)

# 按 tissue 类目 + 病例计数
def cls(t):
    tl = t.lower()
    if "tumour" in tl or "tumor" in tl:
        return "tumour"
    if "normal" in tl:
        return "normal"
    if "adenoma" in tl or "polyp" in tl:
        return "adenoma"
    return "other:" + t
cases = defaultdict(set)
for r in rows:
    cases[cls(r["tissue"])].add(case_key(r["title"]))
print("\n=== 病例级构成（title 首段去重）===", flush=True)
for k, s in cases.items():
    print("  %-14s %d 例: %s" % (k, len(s), sorted(s)), flush=True)

print("\n=== 是否有 adenoma 字样（tissue/title）===", flush=True)
hits = [r for r in rows if re.search(r"adenom|polyp|adenomatous", (r["tissue"] + " " + r["title"]), re.I)]
print("  命中 %d 条" % len(hits))
for r in hits[:20]:
    print("   ", r["gsm"], r["title"], "|", r["tissue"])

json.dump({"rows": rows, "tissue_counts": dict(tc),
           "class_cases": {k: sorted(v) for k, v in cases.items()}},
          open("gse303070_tissue.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved gse303070_tissue.json", flush=True)
