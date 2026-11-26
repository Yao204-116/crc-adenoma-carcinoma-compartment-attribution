# -*- coding: utf-8 -*-
"""补充分析 ⑦：macrodissection 敏感性（B 类遗留 M14b）

审稿人诉求（模拟审稿 M14）：GSE20916 是 macrodissected 的 Carcinoma vs Adenoma 队列——
把整个共识签名建立在含 macro 样本的队列上，签名是否被 macro 混杂驱动？
做法：去掉 GSE20916macro，用其余 3 个发现队列按**完全相同的阈值**
（|meanFC|>1.0 + 全部同向 + 每队列 p<1e-5 + Stouffer q<0.01）重建签名，
对比：① 基因集重叠/Jaccard ② 关键基因（Paneth/分化/间质）留存 ③ 下游三验证队列
无向 AUC 是否保持。

隔离：全部输出写 _loo_sensitivity/{ctrl,no20916}/，绝不覆盖正式 consensus_*.csv。
控制组（4 队列原样重跑）必须逐位复现正式 consensus_signature_tier.csv → 证明改造无害。
用法：python _loo_macro_sensitivity.py   （约 2 分钟）
输出：_loo_macro_sensitivity.json
"""
import csv as _csv
import io, json, os, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
TMP = os.path.join(HERE, "_loo_sensitivity")

# ---------- 0. 搭临时目录 ----------
for sub in ("ctrl", "no20916"):
    os.makedirs(os.path.join(TMP, sub, "02.result"), exist_ok=True)
    for f in ("GSE37364_DEG_Carcinoma_vs_Adenoma.csv",
              "GSE41657_DEG_Carcinoma_vs_Adenoma.csv",
              "GSE20916macro_DEG_Carcinoma_vs_Adenoma.csv"):
        shutil.copy2(os.path.join(HERE, "02.result", f),
                     os.path.join(TMP, sub, "02.result", f))
print("临时目录就绪（02.result 三个 DEG 表已复制）")

src_sig = io.open("consensus_signature.py", encoding="utf-8").read()
src_tier = io.open("consensus_tier_novelty.py", encoding="utf-8").read()
TIER_CUT = src_tier.index("# ---- ② 决定性一问 ----")

MACRO_LINES = '''    ("GSE20916", "Carcinoma vs Adenoma GPL570 (macro)",
     os.path.join(BASE, "02.result/GSE20916macro_DEG_Carcinoma_vs_Adenoma.csv"), "symbol"),
'''


def run_signature(base_dir, drop_macro=False):
    src = src_sig
    if drop_macro:
        # ★ 按行过滤，只删元组两行（docstring 里也有 "GSE20916" 字样，不能误删）
        lines = src.splitlines(keepends=True)
        hit = [i for i, l in enumerate(lines)
               if '("GSE20916"' in l or "GSE20916macro_DEG" in l]
        assert len(hit) == 2, "GSE20916 元组应恰占两行: %s" % hit
        assert hit[1] == hit[0] + 1, "GSE20916 元组两行不连续: %s" % hit
        dropped = "".join(lines[hit[0]:hit[1] + 1])
        assert "GSE41258" not in dropped and "GSE37364" not in dropped
        src = "".join(l for i, l in enumerate(lines) if i not in hit)
    os.environ["AACA_BASE"] = base_dir
    ns = {"__name__": "_cs", "__file__": os.path.join(HERE, "consensus_signature.py")}
    exec(src, ns)


def run_tier(base_dir, three_cohorts=False):
    src = src_tier
    src = src.replace('BASE = r"%s"' % HERE.replace("\\", "/"), 'BASE = r"%s"' % base_dir.replace("\\", "/"))
    if three_cohorts:
        src = src.replace('r["n_cohort"] == 4 and r["agree"] == 4',
                          'r["n_cohort"] == 3 and r["agree"] == 3')
    # ★ 切点必须在替换后的串上重算——BASE 变长会整体移位（曾切进 json.dump 中间）
    cut = src.index("# ---- ② 决定性一问 ----")
    ns = {"__name__": "_tn", "__file__": os.path.join(HERE, "consensus_tier_novelty.py")}
    exec(src[:cut], ns)


def read_tier(path):
    return list(_csv.DictReader(io.open(path, encoding="utf-8-sig")))


# ---------- 1. 控制组（4 队列原样） → 必须复现正式签名 ----------
run_signature(os.path.join(TMP, "ctrl"), drop_macro=False)
run_tier(os.path.join(TMP, "ctrl"), three_cohorts=False)
official = read_tier("consensus_signature_tier.csv")
ctrl = read_tier(os.path.join(TMP, "ctrl", "consensus_signature_tier.csv"))
assert [(r["gene"], r["change"]) for r in official] == \
       [(r["gene"], r["change"]) for r in ctrl], "控制组未复现正式签名！"
print("控制组：4 队列原样重跑 → 逐位复现正式签名（%d 基因）✓" % len(official))

# ---------- 2. LOO：去掉 GSE20916macro ----------
run_signature(os.path.join(TMP, "no20916"), drop_macro=True)
run_tier(os.path.join(TMP, "no20916"), three_cohorts=True)
loo = read_tier(os.path.join(TMP, "no20916", "consensus_signature_tier.csv"))
print("LOO：3 队列（去 macro）同阈值重建 → %d 基因" % len(loo))

# ---------- 3. 对比 ----------
off_set = {r["gene"] for r in official}
loo_set = {r["gene"] for r in loo}
inter = off_set & loo_set
key_genes = ["DEFA5", "DEFA6", "OLFM4", "SPINK1", "LYZ", "REG3G", "PLA2G2A", "ITLN1",
             "MUC2", "TFF3", "FAP", "THBS2", "POSTN", "COL11A1", "INHBA", "ASPN"]
key_status = {g: (g in off_set, g in loo_set) for g in key_genes}
# 每基因 logFC 相关（交集内）
off_map = {r["gene"]: float(r["mean_logFC"]) for r in official}
loo_map = {r["gene"]: float(r["mean_logFC"]) for r in loo}
import statistics
common = sorted(inter)
xs = [off_map[g] for g in common]
ys = [loo_map[g] for g in common]
mx, my = statistics.mean(xs), statistics.mean(ys)
rho = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / (
    (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5)
comp_out = {
    "official_n": len(official), "loo_n": len(loo),
    "overlap": len(inter), "jaccard": round(len(inter) / len(off_set | loo_set), 4),
    "logFC_spearman_pearson_common": round(rho, 4),
    "direction_flip_in_common": sum(1 for g in common if (off_map[g] > 0) != (loo_map[g] > 0)),
    "key_genes": {g: {"official": o, "loo": l} for g, (o, l) in key_status.items()},
}
print("重叠 %d/%d (Jaccard %.3f) | 交集 logFC Pearson r = %.3f | 方向翻转 %d"
      % (len(inter), len(official), comp_out["jaccard"], rho,
         comp_out["direction_flip_in_common"]))
for g, (o, l) in key_status.items():
    if not l:
        print("  ★ 关键基因 LOO 缺失: %s (official=%s)" % (g, o))

# ---------- 4. 下游：两个签名在三验证队列的无向 AUC ----------
_src = io.open("manuscript_facts.py", encoding="utf-8").read()
_ns = {"__name__": "_mfdefs", "__file__": os.path.join(HERE, "manuscript_facts.py")}
exec(_src[:_src.index('print("载入队列")')], _ns)
load_cohort, COHORTS = _ns["load_cohort"], _ns["COHORTS"]
import statistics as st


def score_auc(sig_genes):
    out = {}
    for gse in COHORTS:
        c = load_cohort(gse)
        zmap = c["zmap_adj"]
        gs = [g for g in sig_genes if g in zmap]
        lab = c["labels"]
        sa = [st.mean([zmap[g][i] for g in gs]) for i, l in enumerate(lab) if l == "adenoma"]
        sc = [st.mean([zmap[g][i] for g in gs]) for i, l in enumerate(lab) if l == "carcinoma"]
        out[gse] = round(_ns["auc"](sc, sa), 4)
    return out

auc_official = score_auc(sorted(off_set))
auc_loo = score_auc(sorted(loo_set))
print("下游无向 AUC（癌 vs 腺瘤）：")
print("  官方签名:", auc_official)
print("  LOO 签名:", auc_loo)

json.dump({"meta": {"design": "leave-GSE20916macro-out; identical thresholds but 3/3 cohorts; "
                    "harness validated by control run reproducing the official tier file"},
           "comparison": comp_out,
           "auc_validation": {"official": auc_official, "loo_without_macro": auc_loo}},
          io.open("_loo_macro_sensitivity.json", "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("→ _loo_macro_sensitivity.json")
