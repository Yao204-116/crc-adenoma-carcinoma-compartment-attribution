# -*- coding: utf-8 -*-
"""BMC 稿合规校验：残留引用制式 / 汉字 / 脚手架 / 编号闭环 / 图号顺序 / 词数。"""
import re, json

P = "Manuscript_BMC_Genomics_投稿版.html"
h = open(P, encoding="utf-8").read()
body_html = re.sub(r"<style>.*?</style>", "", h, flags=re.S)
txt = re.sub(r"<[^>]+>", " ", body_html)
for a, b in [("&nbsp;", " "), ("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"),
             ("&mdash;", "—"), ("&middot;", "·"), ("&ndash;", "–")]:
    txt = txt.replace(a, b)
txt = re.sub(r"\s+", " ", txt)

ok = True
def check(name, cond, detail=""):
    global ok
    mark = "PASS" if cond else "FAIL"
    if not cond:
        ok = False
    print(f"[{mark}] {name}" + (f"  {detail}" if detail else ""))

print("=" * 72)
# 1 残留 author-year
ay = re.findall(r"\([^()]{0,60}?(?:19|20)\d{2}[^()]{0,60}?\)", txt)
ay = [x for x in ay if re.search(r"et al\.|[A-Z][a-z]+ & |López", x)]
check("无残留 author-year 引用", not ay, str(ay[:5]))

# 2 残留中文
zh = re.findall(r"[\u4e00-\u9fff]", txt)
check("正文零汉字", not zh, f"{len(zh)} 字" + ("" if not zh else " 例:" + "".join(zh[:20])))

# 3 脚手架
sc = []
for pat in [r"Paragraph \d+ —", r"MS ¶", r"<span class=\"lbl\">", r"¶\d+",
            r"R\d+ ¶", r"必查", r"批注"]:
    m = re.findall(pat, body_html)
    if m:
        sc.append((pat, len(m)))
check("无撰写脚手架标签", not sc, str(sc))

# 4 引用编号闭环（解析 [1, 3-5] 形式的区间）
used = set()
for m in re.finditer(r"\[(\d+(?:\s*[,\-]\s*\d+)*)\]", txt):
    for tok in m.group(1).split(","):
        tok = tok.strip()
        if "-" in tok:
            a, b = [int(x) for x in tok.split("-")]
            used.update(range(a, b + 1))
        elif tok.isdigit():
            used.add(int(tok))
refs = json.load(open("_bmc_references.json", encoding="utf-8"))
nums = {r["num"] for r in refs}
check("正文角标编号 ⊆ 文献表", used <= nums, f"正文用到 {sorted(used)}")
missing = sorted(nums - used)
check("每条文献都被正文引用", not missing, f"未被引用: {missing}")

# 5 图号引用顺序（只看主文 Background → Figure legends 之前）
i_bg = h.find("Background</h2>")
i_fl = h.find("Figure legends</h2>")
body_only = re.sub(r"<[^>]+>", " ", h[i_bg:i_fl])
body_only = re.sub(r"\s+", " ", body_only)
figs = [m.group(0) for m in re.finditer(r"Fig\.\s*S?\d+", body_only)]
seq = []
for f in figs:
    n = re.sub(r"[^A-Za-z0-9]", "", f)
    if n not in seq:
        seq.append(n)
main = [s for s in seq if not s.startswith("FigS")]
check("主图首见序递增", main == sorted(main, key=lambda x: int(re.sub(r"\D", "", x))),
      " → ".join(main))

# 6 摘要
absblk = re.search(r'<div class="ms-para abstract">(.*?)</div>', h, re.S).group(1)
abs_plain = re.sub(r"<[^>]+>", " ", absblk)
for a, b in [("&lt;", "<"), ("&gt;", ">"), ("&amp;", "&"), ("&nbsp;", " ")]:
    abs_plain = abs_plain.replace(a, b)
aw = len(re.findall(r"[A-Za-z][A-Za-z\-']*", abs_plain))
check("摘要 ≤ 350 词（含三段标签）", aw <= 350, f"{aw} 词")
check("摘要三段式标签", all(k in absblk for k in ["Background", "Results", "Conclusions"]))
check("摘要内无文献角标", not re.search(r"\[\d", absblk))

# 7 结构
need = ["Abstract", "Background", "Methods", "Results", "Discussion", "Conclusions",
        "List of abbreviations", "Declarations", "References", "Figure legends",
        "Additional files"]
order = [re.sub(r"<[^>]+>", "", m.group(2)).strip()
         for m in re.finditer(r'<h2 class="chap"><span class="chapno">(\d+)</span>(.*?)</h2>', h)]
check("BMC 节序齐全且顺序正确", order == need, str(order))

# 8 Declarations 8 小节
decl = h[h.find("Declarations"):h.find("References</h2>") if h.find("References</h2>") > h.find("Declarations") else len(h)]
for k in ["Ethics approval and consent to participate", "Consent for publication",
          "Availability of data and materials", "Competing interests", "Funding",
          "Authors' contributions", "Acknowledgements"]:
    if k not in decl:
        check(f"Declarations 含 {k}", False)
check("Declarations 必备小节齐全", all(
    k in decl for k in ["Ethics approval and consent to participate",
                        "Consent for publication", "Availability of data and materials",
                        "Competing interests", "Funding", "Authors' contributions",
                        "Acknowledgements"]))

# 9 图注 ≤ 300 词
idx = [m.start() for m in re.finditer(r'<div class="leg">', h)] + [len(h)]
over = []
for i in range(len(idx) - 1):
    b = h[idx[i]:idx[i + 1]]
    b = b[:b.find("</div>") + 6]
    w = len(re.findall(r"[A-Za-z][A-Za-z\-']*", re.sub(r"<[^>]+>", " ", b)))
    t = re.sub(r"<[^>]+>", "", re.search(r"<h2>(.*?)</h2>", b, re.S).group(1))
    if w > 300:
        over.append((t, w))
check("全部图注 ≤ 300 词", not over, str(over))

# 10 图题 ≤ 15 词
longt = []
for i in range(len(idx) - 1):
    b = h[idx[i]:idx[i + 1]]
    b = b[:b.find("</div>") + 6]
    t = re.sub(r"<[^>]+>", "", re.search(r"<h2>(.*?)</h2>", b, re.S).group(1))
    t = re.sub(r"^Figure\s*S?\d+\.\s*", "", t)
    w = len(re.findall(r"[A-Za-z][A-Za-z\-']*", t))
    if w > 15:
        longt.append((t, w))
check("全部图题 ≤ 15 词", not longt, str(longt))

# 11 关键词
kw = re.search(r'<p class="kw">(.*?)</p>', h, re.S).group(1)
kwn = len([x for x in re.sub(r"<[^>]+>", "", kw).replace("Keywords:", "").split(";") if x.strip()])
check("关键词 3–10 个", 3 <= kwn <= 10, f"{kwn} 个")

# 12 主文词数（Background → Conclusions 结束，不含 Abbreviations/Declarations）
i_bg2 = h.find("Background</h2>")
i_cn = h.find("Conclusions</h2>")
i_ab = h.find("List of abbreviations</h2>")
main_txt = re.sub(r"<[^>]+>", " ", h[i_bg2:i_ab])
mw = len(re.findall(r"[A-Za-z][A-Za-z\-']*", re.sub(r"\s+", " ", main_txt)))
print(f"\n主文（Background→Conclusions）词数: {mw}")
av_txt = re.sub(r"<[^>]+>", " ", h[i_ab:h.find("References</h2>")])
print(f"（Abbreviations + Declarations 另计 {len(re.findall(r'[A-Za-z][A-Za-z-]*', av_txt))} 词）")

# 13 附加文件命名已统一：正文不得再出现 Supplementary Table / Fig. S
resid = {p: len(re.findall(re.escape(p), body_html))
         for p in ["Supplementary Table", "Supplementary File", "Fig. S"]}
check("正文无残留 Supplementary/Fig. S 命名", not any(resid.values()), str(resid))

# 14 附加文件按正文引用顺序编号（1..8 且首见序不回退）
i_af = h.find("Additional files</h2>")
af_body = h[i_bg:i_af] if i_af > i_bg else h
seen = []
for m in re.finditer(r"Additional file (\d+)", af_body):
    n = int(m.group(1))
    if n not in seen:
        seen.append(n)
check("附加文件按正文首见序编号", seen == list(range(1, 9)), f"首见序 {seen}")

# 15 附加文件清单：8 条，且每条都带"文件名/格式/标题/描述"四要素
af_list = h[i_af:h.find("submission version")] if i_af > 0 else ""
items = re.findall(r"<li><b>Additional file (\d+)</b>:(.*?)</li>", af_list, re.S)
lack = [n for n, s in items
        if not re.search(r"\.(xlsx|csv|png|pdf)", s)]
check("附加文件清单 8 条且含文件格式", len(items) == 8 and not lack,
      f"{len(items)} 条" + (f"，缺格式: {lack}" if lack else ""))

# 16 附图图注标题须写成 "Additional file n (Figure Sn)"
hdr = re.findall(r'<span class="fig">Additional file (\d+) \(Figure S(\d)\)\.</span>', h)
check("附图图注标题已挂附加文件号", len(hdr) == 3, str(hdr))

# 17 Availability 段：须给出 10 个 GEO 直链、访问日期，且不得只写"向作者索取"
i_av = h.find("<h3>Availability of data and materials</h3>")
i_av_end = h.find("<h3>Competing interests</h3>")
av = h[i_av:i_av_end] if i_av >= 0 else ""
acc = ["GSE41258", "GSE37364", "GSE20916", "GSE41657", "GSE4183",
       "GSE117606", "GSE15960", "GSE128435", "GSE161277", "GSE303070"]
miss_acc = [a for a in acc if a not in av]
n_url = len(re.findall(r"https://www\.ncbi\.nlm\.nih\.gov/geo/query/acc\.cgi\?acc=GSE\d+", av))
check("Availability 列全 10 个 GEO 永久标识", not miss_acc and n_url == 10,
      f"缺 {miss_acc}；直链 {n_url}/10")
check("Availability 不再依赖“向作者索取”",
      "on reasonable request" not in av and "Code availability" in av)

print("=" * 72)
print("总体：", "全部通过" if ok else "存在未通过项")
