# -*- coding: utf-8 -*-
"""静态检查：出图脚本里是否还残留「上屏用」的中文字符串字面量。

与 `fig_en.finish()` 的运行时守卫互补 —— 两者都过，才算图真的英文化干净：
  * `finish()`（运行时）抓**画布上真出现的**中文：包括 f-string 求值后拼出来的、
    以及从 JSON 数据里带出来的（这两类源码里根本搜不到）。
  * 本脚本（静态）抓**源码里写死的**中文字面量：不用等跑完出图就能发现漏翻。

★ 关键概念：**数据键允许是中文**。
  这些脚本里 `MAT.get("潘氏细胞·4")`、`"基线/结构" in k`、`PAIR["B_paired"]["全签名"]`
  这类**必须**保持中文，翻掉就会静默查空（本仓库踩过：fig5 墨迹覆盖率 47.4% → 7.85%）。
  所以这里用一份**显式白名单**区分「上屏文字」与「数据键」：
  白名单之外的任何中文字面量都报错；以后新增中文字面量时，
  要么翻成英文，要么加进白名单并注明它是数据键。

用法：python check_fig_english.py          # 检查
     python check_fig_english.py -v       # 连白名单命中一起打印
"""
import ast
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SCRIPTS = ["make_figs.py", "make_figs_sc.py", "make_figS1.py", "make_figS3.py"]

CJK = re.compile(r"[\u2e80-\u2eff\u3000-\u303f\u3400-\u4dbf\u4e00-\u9fff"
                 r"\uf900-\ufaff\ufe30-\ufe4f\uff00-\uffef]")

# ---------------------------------------------------------------------------
# 白名单：确认为「数据键」的中文字面量（只按文件 + 字面量内容匹配）
# ---------------------------------------------------------------------------
ALLOWED = {
    "make_figs.py": {
        # fig2(): three_layer_pattern.json 的面板键
        "潘氏细胞·4", "黏液杯状·5", "上皮分化·核心8", "相代谢·5", "纯上皮基线·7", "间质·激活/CAF",
        # fig3(): manuscript_facts.json panel_matrix 的键
        "上调分量(210)", "下调分量(103)", "泛间质·6",
        # fig5(): panel_matrix 的键
        "313签名", "上皮·分化/功能", "间质·基线", "髓系·激活", "Hub基因",
        # fig6(): gse117606_paired_results.json 的键 + 颜色判断用的子串
        "全签名", "上调分量(间质/髓系激活)", "下调分量(上皮功能关闭)", "基线/结构", "阴性",
    },
}


def docstring_nodes(tree):
    """收集所有 docstring 的 Constant 节点 id（这些是文档，允许中文）。"""
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", None)
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                out.add(id(body[0].value))
    return out


def main():
    verbose = "-v" in sys.argv
    total_bad = 0
    for fn in SCRIPTS:
        path = os.path.join(HERE, fn)
        if not os.path.exists(path):
            print("  ?? 缺少脚本 %s" % fn)
            total_bad += 1
            continue
        src = open(path, encoding="utf-8").read()
        tree = ast.parse(src, filename=fn)
        docs = docstring_nodes(tree)
        allowed = ALLOWED.get(fn, set())
        bad, kept = [], []
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
                continue
            if id(node) in docs:
                continue
            if not CJK.search(node.value):
                continue
            if node.value in allowed:
                kept.append(node.value)
            else:
                bad.append((node.lineno, node.value))
        status = "OK " if not bad else "!! "
        print("%s%-18s 残留 %d 处  白名单命中 %d 处" % (status, fn, len(bad), len(kept)))
        for lineno, s in bad:
            print("       line %-4d %r" % (lineno, s))
        if verbose:
            for s in sorted(set(kept)):
                print("       [数据键] %r" % s)
        total_bad += len(bad)
    print()
    if total_bad:
        print("✗ 共有 %d 处上屏中文残留，请翻成英文或加入白名单。" % total_bad)
        return 1
    print("✓ 4 个出图脚本的字符串字面量已无上屏中文（数据键除外）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
