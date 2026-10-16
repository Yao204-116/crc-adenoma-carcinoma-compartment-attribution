# -*- coding: utf-8 -*-
"""GSE303070：FOV 级与供体级组成比较（独立注释，独立平台）。"""
import gzip, csv, json
from collections import defaultdict
import numpy as np
try:
    from scipy.stats import mannwhitneyu
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

P = "sc303070/GSE303070_All_cells_metadata.csv.gz"

EPI = {"Epi"}
STROMA = {"Fibro", "Peri", "SmoothMuscle", "Schwann"}
ENDO = {"Endo"}
IMMUNE = {"Macro", "DC", "Plasma", "Granulo", "Mono", "B", "ILC",
          "TCD8", "TCD4", "Mast", "TZBTB16", "Tgd", "NK"}


def bucket(p):
    if p in EPI:
        return "epi"
    if p in STROMA:
        return "stroma"
    if p in ENDO:
        return "endo"
    if p in IMMUNE:
        return "immune"
    return "other"


fov = defaultdict(lambda: defaultdict(int))
fov_meta = {}
don = defaultdict(lambda: defaultdict(int))
don_meta = {}

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
        b = bucket(row[I["Manual_toplevel_pred"]].strip())
        smp = row[I["Sample"]].strip()
        dv = row[I["Donor"]].strip()
        fk = smp + "|FOV" + row[I["fov"]].strip()
        fov[fk][b] += 1
        fov[fk]["all"] += 1
        fov_meta[fk] = {"donor": dv, "sample": smp, "type": t}
        don[dv]["all"] += 1
        don[dv][b] += 1
        don_meta[dv] = {"sample": smp, "type": t}

print("QC=High 细胞：FOV %d 个，供体 %d 例" % (len(fov), len(don)), flush=True)

def fracs(keys, src):
    out = {"epi": [], "stroma": [], "endo": [], "immune": []}
    for k in keys:
        a = src[k]["all"]
        if a < 20:          # 细胞太少的 FOV 不参与比例估计
            continue
        for b in out:
            out[b].append(src[k].get(b, 0) / a)
    return out

def cmp_level(src, meta, label, do_stats=True):
    print("\n===== %s 级 =====" % label, flush=True)
    tk = [k for k in src if meta[k]["type"] == "T" and src[k]["all"] >= (20 if label == "FOV" else 1)]
    nk = [k for k in src if meta[k]["type"] == "N" and src[k]["all"] >= (20 if label == "FOV" else 1)]
    ft, fn = fracs(tk, src), fracs(nk, src)
    res = {"label": label, "n_T": len(tk), "n_N": len(nk)}
    print("  T: %d 个；N: %d 个" % (len(tk), len(nk)), flush=True)
    print("  %-9s %-22s %-22s %-10s %s" % ("类别", "T 中位[IQR]", "N 中位[IQR]", "Δ(T−N)", "p"), flush=True)
    for b in ["epi", "stroma", "endo", "immune"]:
        a1, a2 = np.array(ft[b]), np.array(fn[b])
        if len(a1) == 0 or len(a2) == 0:
            continue
        m1, m2 = np.median(a1), np.median(a2)
        q1 = np.percentile(a1, [25, 75]); q2 = np.percentile(a2, [25, 75])
        p = float("nan")
        if do_stats and HAVE_SCIPY:
            try:
                p = mannwhitneyu(a1, a2, alternative="two-sided").pvalue
            except Exception:
                pass
        print("  %-9s %.3f [%.3f,%.3f]    %.3f [%.3f,%.3f]    %+.3f    %s"
              % (b, m1, q1[0], q1[1], m2, q2[0], q2[1], m1 - m2,
                 ("%.2e" % p) if p == p else "n/a"), flush=True)
        res[b] = {"T_median": float(m1), "N_median": float(m2), "delta": float(m1 - m2),
                  "T_iqr": [float(q1[0]), float(q1[1])], "N_iqr": [float(q2[0]), float(q2[1])],
                  "p": (float(p) if p == p else None)}
    # 供体级原始占比（不筛选）
    tot = {"T": defaultdict(int), "N": defaultdict(int)}
    for k in src:
        t = meta[k]["type"]
        for b in ["all", "epi", "stroma", "endo", "immune"]:
            tot[t][b] += src[k].get(b, 0)
    res["pooled"] = {t: {b: (tot[t][b] / tot[t]["all"] if tot[t]["all"] else None)
                         for b in ["epi", "stroma", "endo", "immune"]} for t in ("T", "N")}
    print("  合并占比：T epi=%.3f stroma=%.3f endo=%.3f immune=%.3f (n=%d)"
          % (res["pooled"]["T"]["epi"], res["pooled"]["T"]["stroma"], res["pooled"]["T"]["endo"],
             res["pooled"]["T"]["immune"], tot["T"]["all"]), flush=True)
    print("          ：N epi=%.3f stroma=%.3f endo=%.3f immune=%.3f (n=%d)"
          % (res["pooled"]["N"]["epi"], res["pooled"]["N"]["stroma"], res["pooled"]["N"]["endo"],
             res["pooled"]["N"]["immune"], tot["N"]["all"]), flush=True)
    return res

print("HAVE_SCIPY =", HAVE_SCIPY, flush=True)
r_fov = cmp_level(fov, fov_meta, "FOV")
r_don = cmp_level(don, don_meta, "供体")

# 逐供体列出
print("\n===== 逐供体上皮占比（QC=High）=====", flush=True)
for dv in sorted(don, key=lambda d: (don_meta[d]["type"], d)):
    m = don_meta[dv]
    a = don[dv]["all"]
    print("  %-12s %s  细胞 %6d  epi %.3f  stroma %.3f  immune %.3f"
          % (dv, m["type"], a, don[dv]["epi"] / a, don[dv]["stroma"] / a, don[dv]["immune"] / a), flush=True)

json.dump({"fov": r_fov, "donor": r_don,
           "donor_level": {dv: {"type": don_meta[dv]["type"], "all": don[dv]["all"],
                                "epi": don[dv]["epi"] / don[dv]["all"],
                                "stroma": don[dv]["stroma"] / don[dv]["all"],
                                "immune": don[dv]["immune"] / don[dv]["all"]} for dv in don}},
          open("sc303070_composition.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved sc303070_composition.json", flush=True)
