# -*- coding: utf-8 -*-
"""
R3（三层模式 / 腺瘤峰）事实核对：把「双口径筛选」自动化。

背景与动机
----------
本轮 R2 发现了「同一指标两套数字」的坑：three_layer_pattern.json（旧独立脚本
three_layer_check.py 产出，raw 口径、未剔 SSA）与 manuscript_facts.json 的
three_layer（同一管线、含 raw/adj 双口径、已剔 SSA）对同一面板给出不同标签。

R3 的结论必须建立在**双口径都成立**的面板上，所以这里把筛选显式化：
  - 口径 raw  = 队列内 z 化（未做全局漂移校正）
  - 口径 adj  = 队列内 z 化 + 样本级 glev 漂移校正（见 R2 §0）
  - 判据 P/D/A 沿用 manuscript_facts.py 的 pattern()
  - 只有「raw 与 adj 下均为 P」才算稳健腺瘤峰

输出：r3_facts.json（双口径全表 + 筛选结果 + 关键读数）
"""
import json, os, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

F = json.load(open("manuscript_facts.json", encoding="utf-8"))
TL = F["three_layer"]
RAW, ADJ = TL["raw"], TL["adj"]


def pat(no, ad, ca):
    if ad >= no and ad > ca:
        return "P"
    if no > ad > ca:
        return "D"
    return "A"


def label(x):
    return {"P": "P 腺瘤峰", "D": "D 单调降", "A": "A 腺瘤已降"}[x]


print("=" * 104)
print("R3 事实核对：三层模式（正常 → 腺瘤 → 癌）双口径")
print("=" * 104)

names = sorted(set(RAW) | set(ADJ))
rows = {}
print("\n%-22s | %-30s | %-30s | 双口径均P" % ("面板", "原始口径 no/ad/ca → 模式", "校正口径 no/ad/ca → 模式"))
print("-" * 104)
for n in names:
    r = RAW.get(n, {})
    a = ADJ.get(n, {})
    rb = {}
    ab = {}
    for gse in ["GSE4183", "GSE117606", "GSE15960"]:
        if gse in r:
            rb[gse] = pat(r[gse]["normal"], r[gse]["adenoma"], r[gse]["carcinoma"])
        if gse in a:
            ab[gse] = pat(a[gse]["normal"], a[gse]["adenoma"], a[gse]["carcinoma"])
    rs = "".join(rb.get(g, "?") for g in ["GSE4183", "GSE117606", "GSE15960"])
    as_ = "".join(ab.get(g, "?") for g in ["GSE4183", "GSE117606", "GSE15960"])
    ok = (rs == "PPP" and as_ == "PPP")
    rows[n] = {"raw_patterns": rb, "adj_patterns": ab,
               "raw_str": rs, "adj_str": as_, "robust_P_3of3": ok}
    print("%-22s | raw %-6s (4183/117606/15960)      | adj %-6s                    | %s"
          % (n, rs, as_, "★★★ 是" if ok else "—"))

robust = [n for n, v in rows.items() if v["robust_P_3of3"]]
print("\n★ 双口径 + 3/3 队列 均为 P 的面板：%s" % (robust if robust else "（无）"))

# ---------------- 潘氏面板逐队列双口径明细 ----------------
PAN = "潘氏细胞·4"
print("\n" + "=" * 104)
print("★ 潘氏细胞·4（DEFA5 DEFA6 OLFM4 SPINK1）逐队列双口径")
print("=" * 104)
detail = {}
for gse in ["GSE4183", "GSE117606", "GSE15960"]:
    d = {}
    for mode, src in (("raw", RAW), ("adj", ADJ)):
        if gse in src.get(PAN, {}):
            x = src[PAN][gse]
            d[mode] = {"normal": x["normal"], "adenoma": x["adenoma"], "carcinoma": x["carcinoma"],
                       "pattern": pat(x["normal"], x["adenoma"], x["carcinoma"])}
    detail[gse] = d
    if "raw" in d and "adj" in d:
        r0, a0 = d["raw"], d["adj"]
        print("  %-11s raw: no%+7.3f ad%+7.3f ca%+7.3f → %s | adj: no%+7.3f ad%+7.3f ca%+7.3f → %s"
              % (gse, r0["normal"], r0["adenoma"], r0["carcinoma"], label(r0["pattern"]),
                 a0["normal"], a0["adenoma"], a0["carcinoma"], label(a0["pattern"])))

# 腺瘤是否高于癌（三队列共同点，不受"是否回落到正常"影响）
print("\n  腺瘤 vs 癌（峰的核心判据，两口径）：")
for gse in detail:
    for mode in ("raw", "adj"):
        if mode in detail[gse]:
            x = detail[gse][mode]
            print("    %-11s %-4s ad−ca = %+0.3f  (ad %s ca)" % (gse, mode, x["adenoma"] - x["carcinoma"],
                                                                 ">" if x["adenoma"] > x["carcinoma"] else "≤"))

# ---------------- 阴性对照面板 ----------------
print("\n" + "=" * 104)
print("★ 阴性对照：纯上皮基线·7（EPCAM KRT8/18/19 CDH1 CDH17 VIL1）")
print("=" * 104)
for gse in ["GSE4183", "GSE117606", "GSE15960"]:
    for mode, src in (("raw", RAW), ("adj", ADJ)):
        if gse in src.get("纯上皮基线·7", {}):
            x = src["纯上皮基线·7"][gse]
            print("  %-11s %-4s no%+7.3f ad%+7.3f ca%+7.3f → %s"
                  % (gse, mode, x["normal"], x["adenoma"], x["carcinoma"],
                     label(pat(x["normal"], x["adenoma"], x["carcinoma"]))))

# ---------------- 队列设计 ----------------
print("\n" + "=" * 104)
print("★ 队列设计（样本量，供正文与图注）")
print("=" * 104)
for g, v in F["meta"]["cohorts"].items():
    print("  %-11s %-11s %-18s %s" % (g, v["kind"], v["platform"], v["design"]))

out = {"panel_table": rows, "robust_P_panels": robust, "paneth_detail": detail,
       "cohorts": F["meta"]["cohorts"]}
json.dump(out, open("r3_facts.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\n→ r3_facts.json")
