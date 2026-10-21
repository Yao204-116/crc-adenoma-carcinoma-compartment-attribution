# -*- coding: utf-8 -*-
"""
从 GEO 平台 SOFT 文件（*_family.soft.gz）流式提取「探针 → 基因符号」映射表。

为什么需要它：只有 GPL570 / GPL4133 有官方 `.annot.gz`。
  Agilent 的 GPL25373、Illumina 的 GPL21103 等只在 family.soft.gz 里带平台表。
用法：python gpl_soft_map.py GPL25373
输出：geo_cache/GPL25373_map.tsv   （两列：probe_id \t symbol，无表头）
注意：GPL17586 (HTA 2.0) 的 soft 有 7.6 GB，**别碰** —— 遇到需要它的队列直接换。
"""
import gzip
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "cohort_hunt", "geo_cache")

# 列名候选（不同平台大小写/分隔符不一致，实测坑）
SYM_KEYS = ("gene symbol", "gene_symbol", "genesymbol", "symbol", "gene_assignment")


def extract(gpl):
    src = os.path.join(CACHE, "%s_family.soft.gz" % gpl)
    dst = os.path.join(CACHE, "%s_map.tsv" % gpl)
    if not os.path.exists(src) or os.path.getsize(src) < 1000:
        print("%s: 缺少 %s（先跑 geo_dl2.py SOFT:%s）" % (gpl, os.path.basename(src), gpl), flush=True)
        return None
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        print("%s: 已有映射表 %d 行" % (gpl, sum(1 for _ in open(dst, encoding="utf-8"))), flush=True)
        return dst

    hdr = None
    sym_i = id_i = None
    n = 0
    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as f, \
            open(dst, "w", encoding="utf-8") as out:
        in_tbl = False
        for ln in f:
            if ln.startswith("!platform_table_begin"):
                in_tbl = True
                hdr = None
                continue
            if ln.startswith("!platform_table_end"):
                in_tbl = False
                continue
            if not in_tbl:
                continue
            parts = ln.rstrip("\n").split("\t")
            if hdr is None:
                hdr = [p.strip().strip('"').lower() for p in parts]
                for key in SYM_KEYS:
                    if key in hdr:
                        sym_i = hdr.index(key)
                        break
                id_i = hdr.index("id") if "id" in hdr else 0
                if sym_i is None:
                    print("%s: 表头里找不到基因符号列 -> %s" % (gpl, hdr[:12]), flush=True)
                    return None
                print("%s: 表头列 %d，符号列 index=%d (%s)，ID 列 index=%d"
                      % (gpl, len(hdr), sym_i, hdr[sym_i], id_i), flush=True)
                continue
            if len(parts) <= max(sym_i, id_i):
                continue
            sym = parts[sym_i].strip().strip('"')
            if not sym or sym in ("", "---"):
                continue
            pid = parts[id_i].strip().strip('"')
            if not pid:
                continue
            # 多符号用 /// 或 // 分隔，取第一个
            for sep in ("///", "//", ";", ","):
                if sep in sym:
                    sym = sym.split(sep)[0]
                    break
            sym = sym.strip()
            if not sym or sym in ("---", "null"):
                continue
            out.write("%s\t%s\n" % (pid, sym))
            n += 1
    print("%s: 写出 %d 行 -> %s" % (gpl, n, os.path.basename(dst)), flush=True)
    return dst


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["GPL25373"]):
        extract(g)
