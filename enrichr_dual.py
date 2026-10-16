# -*- coding: utf-8 -*-
"""
对「上皮内在层 56」与「间质层 156」成对做 Enrichr 功能富集。

为什么要成对：
  只有把两层放进**同一批数据库**里比，才能证明它们不是同一套生物学。
  若两层富集到的条目完全不重叠 → 「双层结构」被功能富集独立支持。

Enrichr API 的两个坑（本轮踩到）：
  1) addList 必须用 multipart/form-data（urlencoded 一律 400）
  2) 必须先 addList 拿 userListId，再用 userListId 查询（不能直接传基因）

输出：enrichr_results.json（含两层在所有库上的 top 条目）
"""
import json, os, time, urllib.request, urllib.parse, uuid

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "enrichr_results.json")

UP = json.load(open(os.path.join(HERE, "epithelial_core_genes.json"), encoding="utf-8"))
ST = json.load(open(os.path.join(HERE, "stromal_layer_genes.json"), encoding="utf-8"))
print("上皮内在层 %d / 间质层 %d" % (len(UP), len(ST)))

LIBS = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    "Human_Gene_Atlas",
    "PanglaoDB_Augmented_2021",
    "Tabula_Sapiens",
    "Jensen_TISSUES",
]
BASE = "https://maayanlab.cloud/Enrichr"
UA = {"User-Agent": "Mozilla/5.0"}


def mp(fields):
    b = "----wb" + uuid.uuid4().hex
    out = []
    for k, v in fields:
        out.append("--" + b)
        out.append('Content-Disposition: form-data; name="%s"' % k)
        out.append("")
        out.append(v)
    out.append("--" + b + "--")
    out.append("")
    return "\r\n".join(out).encode(), "multipart/form-data; boundary=" + b


def post(path, fields, timeout=60):
    body, ct = mp(fields)
    req = urllib.request.Request(BASE + path, data=body,
                                 headers=dict(UA, **{"Content-Type": ct}))
    for attempt in range(4):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:
            if attempt == 3:
                print("    FAIL %s: %s" % (path, e))
                return None
            time.sleep(2 + 2 * attempt)


def get(path, timeout=60):
    req = urllib.request.Request(BASE + path, headers=UA)
    for attempt in range(4):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=timeout).read())
        except Exception as e:
            if attempt == 3:
                print("    FAIL %s: %s" % (path, e))
                return None
            time.sleep(2 + 2 * attempt)


def add_list(genes, desc):
    return post("/addList", [("list", "\n".join(genes)), ("description", desc)])


def query(uid, lib):
    """Enrichr 返回 [term, p, z, combined, overlapping_genes, adj_p, old_p, old_adjp]"""
    q = urllib.parse.urlencode({"userListId": uid, "backgroundType": lib})
    return get("/enrich?" + q)


RES = {}
for tag, genes in (("上皮内在层", UP), ("间质层", ST)):
    print("\n=== %s（%d 基因）===" % (tag, len(genes)))
    li = add_list(genes, tag)
    if not li:
        continue
    uid = li["userListId"]
    print("  userListId=%s shortId=%s" % (uid, li.get("shortId")))
    RES[tag] = {"n": len(genes), "uid": uid, "libs": {}}
    for lib in LIBS:
        r = query(uid, lib)
        if not r:
            continue
        rows = []
        # ★ 实测格式：[rank, term, p, z, combined, overlapping_genes, adj_p, old_p, old_adjp]
        #   共 9 个字段；早期文档写 8 个会 unpack 报错。
        for row in (r.get(lib) or []):
            row = list(row)
            term, p, z, comb = row[1], row[2], row[3], row[4]
            ov = row[5] if len(row) > 5 else []
            adjp = row[6] if len(row) > 6 else p
            rows.append({"term": term, "p": p, "combined": comb, "adj_p": adjp,
                         "n_overlap": len(ov), "overlap": sorted(ov)[:24]})
        rows.sort(key=lambda x: -(x["combined"] or 0))
        RES[tag]["libs"][lib] = rows[:25]
        print("  ── %s ──" % lib)
        for x in rows[:8]:
            print("    %-62s adjP=%.1e  n=%d  %s"
                  % (x["term"][:62], x["adj_p"] or 1, x["n_overlap"],
                     ",".join(x["overlap"][:8])))
        time.sleep(0.4)

json.dump(RES, open(OUT, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("\n→", OUT)
