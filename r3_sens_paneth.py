# -*- coding: utf-8 -*-
"""
R3 敏感性分析：潘氏面板含干细胞共享基因 OLFM4，必须拆开验证。

问题
----
「潘氏细胞·4」= DEFA5 / DEFA6 / OLFM4 / SPINK1。
其中 OLFM4 同时是肠隐窝基底柱状（干细胞）标志物，其升高可能来自干细胞扩增
而非潘氏细胞。若整面板的腺瘤峰主要由 OLFM4 驱动，则「潘氏程序」的命名不成立。

做法
----
复用 manuscript_facts.py 的队列加载（exec 前 227 行，避开其模块级产出覆盖），
对以下子面板重算三层模式（正常/腺瘤/癌）双口径：
  - 全 4 基因（原面板）
  - DEFA5 + DEFA6（纯潘氏）
  - OLFM4（干细胞/潘氏共享）
  - SPINK1（潘氏分泌）
  - DEFA5 + DEFA6 + SPINK1（去 OLFM4）

判据与 R3 主表一致：P=腺瘤峰(ad>=no 且 ad>ca)、D=单调降、A=其余。
"""
import os, statistics, json

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

src = open("manuscript_facts.py", encoding="utf-8").read().split("\n")
g = {"__name__": "mf_head", "__file__": os.path.join(HERE, "manuscript_facts.py")}
exec("\n".join(src[:227]), g)          # 到 IDX 构建完成为止
C, IDX = g["C"], g["IDX"]


def pat(no, ad, ca):
    if ad >= no and ad > ca:
        return "P 腺瘤峰"
    if no > ad > ca:
        return "D 单调降"
    return "A 腺瘤已降"


SUBS = {
    "全4基因(DEFA5/6+OLFM4+SPINK1)": ["DEFA5", "DEFA6", "OLFM4", "SPINK1"],
    "纯潘氏(DEFA5+DEFA6)": ["DEFA5", "DEFA6"],
    "仅OLFM4(干细胞共享)": ["OLFM4"],
    "仅SPINK1": ["SPINK1"],
    "去OLFM4(DEFA5/6+SPINK1)": ["DEFA5", "DEFA6", "SPINK1"],
}

print("=" * 104)
print("R3 敏感性：潘氏面板拆解（三层模式，双口径）")
print("=" * 104)
out = {}
for name, genes in SUBS.items():
    print("\n── %s ──" % name)
    out[name] = {}
    for mode in ("raw", "adj"):
        line = "   %-4s" % mode
        for gse in ["GSE4183", "GSE117606", "GSE15960"]:
            zm = C[gse]["zmap" if mode == "raw" else "zmap_adj"]
            gs = [x for x in genes if x in zm]
            if not gs:
                line += "  %s: 无" % gse
                continue
            med = {}
            for grp in ("no", "ad", "ca"):
                vals = []
                for i in IDX[gse][grp]:
                    zz = [zm[x][i] for x in gs if i < len(zm[x])]
                    if zz:
                        vals.append(statistics.mean(zz))
                med[grp] = statistics.median(vals)
            p = pat(med["no"], med["ad"], med["ca"])
            out[name].setdefault(gse, {})[mode] = {
                "normal": med["no"], "adenoma": med["ad"], "carcinoma": med["ca"],
                "pattern": p, "n_genes": len(gs)}
            line += "  %s[%d基因] %+6.3f→%+6.3f→%+6.3f %s" % (
                gse.replace("GSE", ""), len(gs), med["no"], med["ad"], med["ca"], p)
        print(line)

# 汇总：双口径是否均 P、且是否满足 no<ca<ad
print("\n" + "=" * 104)
print("汇总判读")
print("=" * 104)
for name, d in out.items():
    okP = all(d[gse][m]["pattern"] == "P 腺瘤峰" for gse in d for m in d[gse])
    incr = all(d[gse][m]["normal"] < d[gse][m]["carcinoma"] < d[gse][m]["adenoma"]
               for gse in d for m in d[gse])
    print("  %-30s 双口径+3队列全P: %-4s   全局 no<ca<ad: %s" % (name, "是" if okP else "否", "是" if incr else "否"))

json.dump(out, open("r3_sens_paneth.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\n→ r3_sens_paneth.json")
