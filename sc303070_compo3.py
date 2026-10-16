# -*- coding: utf-8 -*-
"""GSE303070 组成分析·终版（修正上一版 bug：donor 键必须带 sample_type）。
键统一为 (donor, sample_type, sample)。"""
import gzip, csv, json
from collections import defaultdict
import numpy as np
from scipy.stats import mannwhitneyu, wilcoxon

P = "sc303070/GSE303070_All_cells_metadata.csv.gz"
EPI = {"Epi"}
STROMA = {"Fibro", "Peri", "SmoothMuscle", "Schwann"}
ENDO = {"Endo"}
IMMUNE = {"Macro", "DC", "Plasma", "Granulo", "Mono", "B", "ILC",
          "TCD8", "TCD4", "Mast", "TZBTB16", "Tgd", "NK"}


def bucket(p):
    return ("epi" if p in EPI else "stroma" if p in STROMA else
            "endo" if p in ENDO else "immune" if p in IMMUNE else "other")


cell = defaultdict(lambda: defaultdict(int))   # (donor,type,sample) -> bucket -> n
fov = defaultdict(lambda: defaultdict(int))    # (donor,type,sample,fov) -> bucket -> n
fmeta = {}

with gzip.open(P, "rt", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f)
    cols = next(rd)
    I = {c: i for i, c in enumerate(cols)}
    for row in rd:
        if len(row) < len(cols):
            continue
        if row[I["Quality"]].strip().lower() != "high":
            continue
        st = row[I["Sample_type"]].strip().upper()
        t = "T" if st.startswith("T") else ("N" if st.startswith("N") else "?")
        if t == "?":
            continue
        d = row[I["Donor"]].strip()
        s = row[I["Sample"]].strip()
        b = bucket(row[I["Manual_toplevel_pred"]].strip())
        k = (d, t, s)
        cell[k]["all"] += 1
        cell[k][b] += 1
        fk = k + (row[I["fov"]].strip(),)
        fov[fk]["all"] += 1
        fov[fk][b] += 1

OUT = {"per_sample": {}, "paired": {}, "fov_global": {}, "fov_niche_confound": {}}

print("=== 逐样本组成（QC=High）===", flush=True)
tbl = []
for k in sorted(cell, key=lambda x: (x[1], x[0])):
    a = cell[k]["all"]
    if a == 0:
        continue
    fr = {b: cell[k][b] / a for b in ("epi", "stroma", "endo", "immune")}
    tbl.append((k, a, fr))
    OUT["per_sample"]["%s|%s" % (k[0], k[1])] = {"sample": k[2], "n": a, **fr}
    print("  %-10s %s %-7s n=%6d epi %.3f stroma %.3f endo %.3f immune %.3f"
          % (k[0], k[1], k[2], a, fr["epi"], fr["stroma"], fr["endo"], fr["immune"]), flush=True)

# 配对供体（同一 donor 同时有 T 与 N 样本）
byt = defaultdict(dict)
for k, a, fr in tbl:
    byt[k[0]][k[1]] = (k[2], a, fr)
paired = {d: v for d, v in byt.items() if "T" in v and "N" in v}
print("\n=== 配对供体（T 与 N 同源）n=%d ===" % len(paired), flush=True)
if paired:
    print("  %-10s %-28s %-28s %s" % ("供体", "T: n(epi%)", "N: n(epi%)", "Δepi(T−N)"), flush=True)
    d_epi, d_str, d_imm = [], [], []
    for d, v in sorted(paired.items()):
        sT, nT, fT = v["T"]; sN, nN, fN = v["N"]
        d_epi.append(fT["epi"] - fN["epi"]); d_str.append(fT["stroma"] - fN["stroma"]); d_imm.append(fT["immune"] - fN["immune"])
        print("  %-10s %-28s %-28s %+.3f" % (d, "%s: %d (%.3f)" % (sT, nT, fT["epi"]),
                                             "%s: %d (%.3f)" % (sN, nN, fN["epi"]), fT["epi"] - fN["epi"]), flush=True)
    OUT["paired"] = {"n": len(paired), "d_epi": d_epi, "d_stroma": d_str, "d_immune": d_imm,
                     "median_d_epi": float(np.median(d_epi))}
    print("  配对 Δepi 中位 %+.3f（同向 %d/%d）" % (np.median(d_epi), sum(1 for x in d_epi if x > 0), len(d_epi)), flush=True)

# FOV 级（细胞≥20），全样本池化
FT, FN = [], []
for fk in fov:
    if fov[fk]["all"] < 20:
        continue
    d, t, s, fv = fk
    fr = {b: fov[fk][b] / fov[fk]["all"] for b in ("epi", "stroma", "endo", "immune")}
    (FT if t == "T" else FN).append(fr)
print("\n=== FOV 级（≥20 细胞）T n=%d / N n=%d ===" % (len(FT), len(FN)), flush=True)
for b in ("epi", "stroma", "endo", "immune"):
    a1 = np.array([x[b] for x in FT]); a2 = np.array([x[b] for x in FN])
    p = mannwhitneyu(a1, a2, alternative="two-sided").pvalue
    print("  %-8s T %.3f [%.3f,%.3f]  N %.3f [%.3f,%.3f]  Δ%+.3f  p=%.2e"
          % (b, np.median(a1), *np.percentile(a1, [25, 75]), np.median(a2),
             *np.percentile(a2, [25, 75]), np.median(a1) - np.median(a2), p), flush=True)
    OUT["fov_global"][b] = {"T_median": float(np.median(a1)), "N_median": float(np.median(a2)),
                            "delta": float(np.median(a1) - np.median(a2)), "p": float(p),
                            "n_T": len(a1), "n_N": len(a2)}

json.dump(OUT, open("sc303070_composition_final.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\nsaved sc303070_composition_final.json", flush=True)
