# -*- coding: utf-8 -*-
"""诊断：GSE303070 的 FOV 是否选择性取样（niche 分布 × 样本类型），
以及 Sample × Donor × Sample_type 的对应关系。"""
import gzip, csv, json
from collections import Counter, defaultdict

P = "sc303070/GSE303070_All_cells_metadata.csv.gz"
niche = defaultdict(Counter)       # sample_type -> niche -> cells
so = defaultdict(Counter)          # sample -> donor
stype = {}                         # sample -> type
fov_niche = defaultdict(Counter)   # sample_type -> niche -> #FOV
seen = set()

with gzip.open(P, "rt", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f)
    cols = next(rd)
    I = {c: i for i, c in enumerate(cols)}
    for row in rd:
        if len(row) < len(cols):
            continue
        s = row[I["Sample"]].strip()
        d = row[I["Donor"]].strip()
        st = row[I["Sample_type"]].strip().upper()
        t = "T" if st.startswith("T") else ("N" if st.startswith("N") else st)
        nz = row[I["niches"]].strip()
        fv = row[I["fov"]].strip()
        stype[s] = t
        so[s][d] += 1
        niche[t][nz] += 1
        k = (s, fv)
        if k not in seen:
            seen.add(k)
            fov_niche[t][nz] += 1

print("=== Sample → Donor → 类型 ===", flush=True)
for s in sorted(so, key=lambda x: (stype[x], x)):
    print("  %-12s %-9s 细胞 %6d  donors=%s" % (s, stype[s], sum(so[s].values()),
          ", ".join("%s(%d)" % (k, v) for k, v in so[s].most_common())), flush=True)

print("\n=== Sample_type × niches（细胞数）===", flush=True)
for t in ("T", "N"):
    tot = sum(niche[t].values())
    print(" %s 合计 %d" % (t, tot), flush=True)
    for k, v in niche[t].most_common(20):
        print("    %-26s %7d  %5.1f%%" % (k, v, 100 * v / tot), flush=True)

print("\n=== Sample_type × niches（FOV 数）===", flush=True)
for t in ("T", "N"):
    tot = sum(fov_niche[t].values())
    print(" %s 合计 %d FOV" % (t, tot), flush=True)
    for k, v in fov_niche[t].most_common(20):
        print("    %-26s %4d FOV  %5.1f%%" % (k, v, 100 * v / tot), flush=True)

json.dump({"sample_donor": {s: dict(so[s]) for s in so}, "stype": stype,
           "cells_by_niche": {t: dict(niche[t]) for t in niche},
           "fov_by_niche": {t: dict(fov_niche[t]) for t in fov_niche}},
          open("sc303070_niche.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved sc303070_niche.json", flush=True)
