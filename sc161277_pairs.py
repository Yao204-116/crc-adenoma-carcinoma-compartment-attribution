# -*- coding: utf-8 -*-
"""GSE161277 配对分解正式分析 v2：腺瘤期分化/潘氏信号 = 组成还是状态？

v2 变更（2026-10-15）：
  - 同一 (patient, stage) 有多个样本时（Patient3 有两个腺瘤样本 GSM4904242/243）
    不再被 dict 覆盖，改为**按细胞数加权合并**：
      T 层（全细胞均值）权重 = n_cells_qc；U 层（上皮内均值）权重 = n_epi
    （每细胞 CP10K+log1p 归一化后，跨细胞均值可精确加权合并）
  - 输出保留 per-sample 行（suppl/绘图用）+ 合并行（分解用）
  - 打印方向一致性统计
"""
import json, math, os

BASE = os.path.dirname(os.path.abspath(__file__))
rs = json.load(open(os.path.join(BASE, "sc161277_scan.json"), encoding="utf-8"))

# ---------- 1. 按 (patient, stage) 加权合并 ----------
SIG_NAMES = ["stem", "prolif", "diff_abs", "goblet", "paneth"]
U_TOP = ["immune_U_epi", "stroma_U_epi", "epi_core_U_epi"]
T_TOP = ["immune_T_all", "stroma_T_all", "epi_core_T_all"]

groups = {}
for r in rs:
    groups.setdefault((r["patient"], r["stage"]), []).append(r)


def wavg(rows, key, wkey):
    num = den = 0.0
    for r in rows:
        w = r[wkey]
        num += r[key] * w
        den += w
    return num / den if den else 0.0


def make_row(rows):
    """把一个 (patient, stage) 下的所有样本合并成一行"""
    n_qc = sum(r["n_cells_qc"] for r in rows)
    n_epi = sum(r["_epi_idx_n"] for r in rows)
    m = {"patient": rows[0]["patient"], "stage": rows[0]["stage"],
         "n_cells_raw": sum(r["n_cells_raw"] for r in rows),
         "n_cells_qc": n_qc, "_epi_idx_n": n_epi,
         "pct_epi_of_qc": n_epi / n_qc if n_qc else 0.0,
         "n_samples": len(rows), "gsms": [r["gsm"] for r in rows]}
    for name in SIG_NAMES:
        m[name] = {
            "T_all_cells": wavg([{**r[name], "w": r["n_cells_qc"]} for r in rows], "T_all_cells", "w"),
            "U_epi_cells": wavg([{**r[name], "w": r["_epi_idx_n"]} for r in rows], "U_epi_cells", "w"),
            "T_non_epi": wavg([{**r[name], "w": r["n_cells_qc"]} for r in rows], "T_non_epi", "w"),
            "frac_epi_pos": wavg([{**r[name], "w": r["_epi_idx_n"]} for r in rows], "frac_epi_pos", "w"),
        }
    for k in U_TOP:
        m[k] = wavg([{**r, "w": r["_epi_idx_n"]} for r in rows], k, "w")
    for k in T_TOP:
        m[k] = wavg([{**r, "w": r["n_cells_qc"]} for r in rows], k, "w")
    m["genes_U_epi"] = {g: wavg([{**r["genes_U_epi"], "w": r["_epi_idx_n"]} for r in rows], g, "w")
                        for g in rows[0]["genes_U_epi"]}
    return m


S = {k: make_row(v) for k, v in groups.items()}
print("=" * 100)
print("合并结果（同一 patient/stage 多样本加权合并）")
print("=" * 100)
for (p, s), v in sorted(S.items()):
    if v["n_samples"] > 1:
        print("  %-10s %-12s  %d 样本合并  n_qc=%d n_epi=%d epi%%=%.1f  gsm=%s"
              % (p, s, v["n_samples"], v["n_cells_qc"], v["_epi_idx_n"],
                 v["pct_epi_of_qc"] * 100, ",".join(v["gsms"])))

# 纯潘氏（去 LYZ）
PANETH_CLEAN = ["DEFA5", "DEFA6", "REG3G", "PLA2G2A", "ITLN1"]


def compose_clean_paneth(r):
    pn = r["paneth"]["U_epi_cells"] * 6
    d5, d6, lyz = r["genes_U_epi"]["DEFA5"], r["genes_U_epi"]["DEFA6"], r["genes_U_epi"]["LYZ"]
    rest = (pn - d5 - d6 - lyz) / 3.0
    return (d5 + d6 + rest) / 3.0


def lg(a, b):
    return math.log(a / b) if a > 0 and b > 0 else None


out = {"decomp": [], "defa": [], "defa_persample": [], "notes": []}

print()
print("=" * 100)
print("配对分解：Δlog T(组织层,≈LCM) = Δlog P(上皮占比) + Δlog U(上皮内状态)   [合并行]")
print("=" * 100)
print("%-10s %-16s %8s %8s %8s | %8s %8s %8s %8s %8s"
      % ("patient", "transition", "dlogT", "dlogP", "dlogU", "P_a%", "P_b%", "U_a", "U_b", "n_epi(b)"))
for p in ["Patient0", "Patient1", "Patient2", "Patient3"]:
    have_ad = (p, "adenoma") in S
    have_ca = (p, "carcinoma") in S
    have_no = (p, "normal") in S
    pairs = []
    if have_ad and have_ca:
        pairs.append(("adenoma", "carcinoma"))
    if have_no and have_ad:
        pairs.append(("normal", "adenoma"))
    for a, b in pairs:
        ra, rb = S[(p, a)], S[(p, b)]
        Pa, Pb = ra["pct_epi_of_qc"], rb["pct_epi_of_qc"]
        for sig in ["paneth", "diff_abs", "goblet", "stem"]:
            Ta, Tb = ra[sig]["T_all_cells"], rb[sig]["T_all_cells"]
            Ua, Ub = ra[sig]["U_epi_cells"], rb[sig]["U_epi_cells"]
            dt, dp, du = lg(Tb, Ta), lg(Pb, Pa), lg(Ub, Ua)
            if dt is None:
                continue
            out["decomp"].append({"patient": p, "a": a, "b": b, "sig": sig,
                                  "T_a": Ta, "T_b": Tb, "U_a": Ua, "U_b": Ub,
                                  "P_a": Pa, "P_b": Pb, "dlogT": dt, "dlogP": dp, "dlogU": du,
                                  "n_epi_a": ra["_epi_idx_n"], "n_epi_b": rb["_epi_idx_n"],
                                  "n_samples_a": ra["n_samples"], "n_samples_b": rb["n_samples"]})
            if sig == "paneth":
                print("%-10s %-16s %8.2f %8.2f %8.2f | %7.1f%% %7.1f%% %8.3f %8.3f %8d"
                      % (p, a + "->" + b[:4], dt, dp or 0, du or 0, Pa * 100, Pb * 100, Ua, Ub, rb["_epi_idx_n"]))

print()
print("=" * 100)
print("DEFA5 / DEFA6 上皮内表达（U 层，合并行）")
print("=" * 100)
print("%-18s %8s %8s %8s %8s %8s %6s" % ("sample", "DEFA5", "DEFA6", "panethU", "panU_cln", "epi_frac", "n_epi"))
for p in ["Patient0", "Patient1", "Patient2", "Patient3"]:
    for s in ["normal", "para-cancer", "adenoma", "carcinoma"]:
        if (p, s) not in S:
            continue
        r = S[(p, s)]
        ge = r["genes_U_epi"]
        pc = compose_clean_paneth(r)
        print("%-18s %8.3f %8.3f %8.3f %8.3f %6.1f%% %6d" %
              (p + " " + s, ge["DEFA5"], ge["DEFA6"], r["paneth"]["U_epi_cells"], pc,
               r["pct_epi_of_qc"] * 100, r["_epi_idx_n"]))
        out["defa"].append({"patient": p, "stage": s, "DEFA5_U": ge["DEFA5"], "DEFA6_U": ge["DEFA6"],
                            "paneth_U": r["paneth"]["U_epi_cells"], "paneth_clean_U": pc,
                            "epi_frac": r["pct_epi_of_qc"], "n_qc": r["n_cells_qc"],
                            "n_epi": r["_epi_idx_n"], "n_samples": r["n_samples"], "merged": r["n_samples"] > 1})

# per-sample（供图 8 与附表）
for r0 in rs:
    ge = r0["genes_U_epi"]
    out["defa_persample"].append({"patient": r0["patient"], "stage": r0["stage"], "gsm": r0["gsm"],
                                  "DEFA5_U": ge["DEFA5"], "DEFA6_U": ge["DEFA6"],
                                  "paneth_U": r0["paneth"]["U_epi_cells"],
                                  "paneth_clean_U": compose_clean_paneth(r0),
                                  "epi_frac": r0["pct_epi_of_qc"], "n_epi": r0["_epi_idx_n"]})

# ---------- 方向一致性统计 ----------
print()
print("=" * 100)
print("方向一致性（DEFA6 / 去 LYZ 潘氏 U 层）")
print("=" * 100)
def cmp_stat(field, cmp_pairs, label):
    for a, b, ttl in cmp_pairs:
        hit = tot = 0
        det = []
        for p in ["Patient1", "Patient2", "Patient3"]:
            if (p, a) in S and (p, b) in S:
                va, vb = S[(p, a)][field], S[(p, b)][field]
                tot += 1
                ok = va > vb
                hit += ok
                det.append("%s %.3f%s%.3f" % (p[-1], va, ">" if ok else "<", vb))
        print("  %-26s %s  %d/%d   %s" % (ttl, label, hit, tot, "; ".join(det)))

def g6(r):
    return r["genes_U_epi"]["DEFA6"]

def clean(r):
    return compose_clean_paneth(r)

pairs_cmp = [("adenoma", "normal", "腺瘤 > 正常"), ("adenoma", "carcinoma", "腺瘤 > 癌")]
for fld, fname in [(None, "DEFA6"), (None, "去LYZ潘氏")]:
    pass
# DEFA6
for a, b, ttl in pairs_cmp:
    hit = tot = 0; det = []
    for p in ["Patient1", "Patient2", "Patient3"]:
        if (p, a) in S and (p, b) in S:
            va, vb = g6(S[(p, a)]), g6(S[(p, b)]); tot += 1
            ok = va > vb; hit += ok
            det.append("%s %.3f%s%.3f" % (p[-1], va, ">" if ok else "<", vb))
    print("  %-26s DEFA6      %d/%d   %s" % (ttl, hit, tot, "; ".join(det)))
for a, b, ttl in pairs_cmp:
    hit = tot = 0; det = []
    for p in ["Patient1", "Patient2", "Patient3"]:
        if (p, a) in S and (p, b) in S:
            va, vb = clean(S[(p, a)]), clean(S[(p, b)]); tot += 1
            ok = va > vb; hit += ok
            det.append("%s %.3f%s%.3f" % (p[-1], va, ">" if ok else "<", vb))
    print("  %-26s 去LYZ潘氏  %d/%d   %s" % (ttl, hit, tot, "; ".join(det)))

# stem U 层 normal->adenoma
hit = tot = 0; det = []
for p in ["Patient1", "Patient2", "Patient3"]:
    if (p, "normal") in S and (p, "adenoma") in S:
        va = S[(p, "normal")]["stem"]["U_epi_cells"]; vb = S[(p, "adenoma")]["stem"]["U_epi_cells"]
        hit += vb > va; tot += 1
        det.append("%s %+.3f" % (p[-1], vb - va))
print("  %-26s stem U    %d/%d   %s" % ("normal→腺瘤 上升", hit, tot, "; ".join(det)))

json.dump(out, open(os.path.join(BASE, "sc161277_pairs.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
print("\nsaved sc161277_pairs.json")
