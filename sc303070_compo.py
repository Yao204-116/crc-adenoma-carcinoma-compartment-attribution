# -*- coding: utf-8 -*-
"""GSE303070 空间层组成交叉验证：
用作者自己的细胞类型注释（Manual_toplevel_pred）比较肿瘤 vs 正常的
上皮占比 / 间质占比 / 免疫占比——独立于我们的 LCM 与签名体系。
单位：① 供体级（每例一个值，取该例所有 QC 细胞合并）② FOV 级（352 视野）。
"""
import gzip, csv, json, sys
from collections import Counter, defaultdict

P = "sc303070/GSE303070_All_cells_metadata.csv.gz"

idx = {}
cnt_pred = Counter()
cnt_stype = Counter()
cnt_quality = Counter()
# 供体级累加
don = defaultdict(lambda: defaultdict(int))
don_meta = {}
# FOV 级
fov = defaultdict(lambda: defaultdict(int))
fov_meta = {}

with gzip.open(P, "rt", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f)
    cols = next(rd)
    for i, c in enumerate(cols):
        idx[c] = i
    I = idx
    n = 0
    for row in rd:
        n += 1
        if len(row) < len(cols):
            continue
        st = row[I["Sample_type"]].strip()
        q = row[I["Quality"]].strip()
        pred = row[I["Manual_toplevel_pred"]].strip()
        smp = row[I["Sample"]].strip()
        dv = row[I["Donor"]].strip()
        fv = row[I["fov"]].strip()
        cnt_stype[st] += 1
        cnt_quality[q] += 1
        cnt_pred[pred] += 1
        # 只统计高质量细胞
        if q.lower() != "high":
            continue
        key = "T" if st.upper().startswith("T") else ("N" if st.upper().startswith("N") else st)
        don[dv][key + "|all"] += 1
        don[dv][key + "|" + pred] += 1
        don_meta[dv] = {"sample": smp, "type": key}
        fk = smp + "|" + fv
        fov[fk][pred] += 1
        fov[fk]["__all"] += 1
        fov_meta[fk] = {"donor": dv, "sample": smp, "type": key}

print("总细胞行数:", n, flush=True)
print("\n=== Sample_type 分布 ===", flush=True)
for k, v in cnt_stype.most_common():
    print("  %-12s %d" % (k, v), flush=True)
print("\n=== Quality 分布 ===", flush=True)
for k, v in cnt_quality.most_common():
    print("  %-12s %d" % (k, v), flush=True)
print("\n=== Manual_toplevel_pred 分布（全部细胞）===", flush=True)
for k, v in cnt_pred.most_common(30):
    print("  %-28s %d" % (k, v), flush=True)

json.dump({"n_rows": n, "sample_type": dict(cnt_stype), "quality": dict(cnt_quality),
           "pred": dict(cnt_pred)}, open("sc303070_meta_counts.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\nsaved sc303070_meta_counts.json", flush=True)
