# -*- coding: utf-8 -*-
"""生成《BMC Genomics 投稿准备报告》HTML。"""
import re, json, os, html as H

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)

MS = "Manuscript_BMC_Genomics_投稿版.html"
OUT = "BMC投稿准备报告_BMC_Genomics.html"

refs = json.load(open("_bmc_references.json", encoding="utf-8"))
refs = sorted(refs, key=lambda r: r["num"])
ms = open(MS, encoding="utf-8").read()

# 每条文献在正文中的落点章节
CHAP_RE = re.compile(r'<h2 class="chap"><span class="chapno">(\d+)</span>(.*?)</h2>')
ms_marks = list(CHAP_RE.finditer(ms))
sections = []
for i, m in enumerate(ms_marks):
    end = ms_marks[i + 1].start() if i + 1 < len(ms_marks) else len(ms)
    sections.append((re.sub(r"<[^>]+>", "", m.group(2)).strip(), ms[m.end():end]))

def where(n):
    """找出编号 n 出现在哪些章节（解析 [a, b-c] 形式的角标）"""
    hits = []
    for name, body in sections:
        if name in ("References", "Figure legends", "Additional files"):
            continue
        txt = re.sub(r"<[^>]+>", " ", body)
        found = False
        for m in re.finditer(r"\[(\d[\d,\s\-]*)\]", txt):
            for tok in m.group(1).split(","):
                tok = tok.strip()
                if "-" in tok:
                    a, b = [int(x) for x in tok.split("-")]
                    if a <= n <= b:
                        found = True
                elif tok.isdigit() and int(tok) == n:
                    found = True
        if found:
            hits.append(name)
    return "、".join(hits) if hits else "—"

# 摘要词数
absblk = re.search(r'<div class="ms-para abstract">(.*?)</div>', ms, re.S).group(1)
absp = re.sub(r"<[^>]+>", " ", absblk)
for a, b in [("&lt;", "<"), ("&gt;", ">")]:
    absp = absp.replace(a, b)
ABS_WORDS = len(re.findall(r"[A-Za-z][A-Za-z\-']*", absp))

# 主文词数（Background → Conclusions，不含 Abbreviations / Declarations）
i_bg = ms.find("Background</h2>")
i_ref = ms.find("List of abbreviations</h2>")
MAIN_WORDS = len(re.findall(r"[A-Za-z][A-Za-z\-']*", re.sub(r"<[^>]+>", " ", ms[i_bg:i_ref])))

FIG_MAP = [("Fig. 1", "Study design and analytical framework", "Background", "Fig. 1", "不变"),
           ("Fig. 6", "56 matched patient pairs", "Results（R1 ¶3）", "Fig. 2", "改号 6→2"),
           ("Fig. 3", "2×2 decomposition", "Results", "Fig. 3", "不变"),
           ("Fig. 4", "Per-gene attribution scatter", "Results", "Fig. 4", "不变"),
           ("Fig. 5", "Panel × cohort discrimination", "Results", "Fig. 5", "不变"),
           ("Fig. 2", "Three-layer pattern", "Results（R3）", "Fig. 6", "改号 2→6"),
           ("Fig. 7", "Single-cell decomposition", "Results（R4）", "Fig. 7", "不变"),
           ("Fig. 8", "DEFA5/DEFA6 in epithelial cells", "Results（R4）", "Fig. 8", "不变"),
           ("Fig. S1", "Global drift diagnostics", "Results", "Fig. S1", "不变"),
           ("Fig. S2", "Dual-caliber screening matrix", "Results", "Fig. S2", "不变"),
           ("Fig. S3", "GSE303070 不可作基准", "Results（R5）", "Fig. S3", "不变")]

REF_ROWS = "".join(
    f'<tr><td class="c"><b>{r["num"]}</b></td><td class="cit">{H.escape(r["vancouver"])}</td>'
    f'<td class="c">{r["pmid"] if not r.get("preprint") else "preprint"}</td>'
    f'<td>{ ("原始" if r["num"]>=17 else "文献") }</td>'
    f'<td>{where(r["num"])}</td></tr>' for r in refs)

FIG_ROWS = "".join(
    f'<tr><td>{a}</td><td>{b}</td><td>{c}</td>'
    f'<td class="c"><b>{d}</b></td><td class="c">{e}</td></tr>' for a, b, c, d, e in FIG_MAP)

# ---- 附加文件映射：直接从稿件解析，保证与成品一致 ----
i_af = ms.find("Additional files</h2>")
AF_ITEMS = []
for m in re.finditer(
        r"<li><b>Additional file (\d+)</b>: (.*?) \((.*?)\) &mdash; (.*?)</li>", ms[i_af:]):
    AF_ITEMS.append({"num": int(m.group(1)), "label": m.group(2),
                     "fmt": m.group(3), "title": m.group(4)})


def af_where(n):
    """该附加文件在正文/图注中的落点"""
    hits = []
    for name, body in sections:
        if name in ("References", "Additional files"):
            continue
        if re.search(rf"Additional file {n}\b", re.sub(r"<[^>]+>", " ", body)):
            hits.append("图注" if name == "Figure legends" else name)
    return "、".join(hits) if hits else "—"


# 上传文件名建议：表格合并为单个 xlsx，图用 png
def af_upload(a):
    ext = "png" if a["label"].startswith("Figure") else "xlsx"
    return f"Additional_file_{a['num']}.{ext}"


# 原内部编号（人工对照，用于说明"AF 号 ≠ S 号"的原因）
AF_SRC = {1: "suppl/S1_313gene_attribution.csv", 2: "suppl/S2_panel_gene_membership.csv",
          3: "suppl/S4a/S4b_sc161277_*.csv", 4: "suppl/S5_four_cohort_performance.csv",
          5: "figs/figS1_drift_diagnostics.png",
          6: "figs/figS2_double_caliber_screening_matrix.png",
          7: "suppl/S3_panel_pattern_matrix.csv + S3b_paneth_sensitivity.csv",
          8: "figs/figS3_GSE303070_sampling_bias.png"}

AF_ROWS = "".join(
    f'<tr><td class="c"><b>Additional file {a["num"]}</b></td>'
    f'<td class="c">{a["label"]}</td>'
    f'<td>{a["title"]}</td>'
    f'<td class="c">{a["fmt"]}</td>'
    f'<td class="c">{af_where(a["num"])}</td>'
    f'<td class="mono">{af_upload(a)}</td>'
    f'<td class="mono">{AF_SRC[a["num"]]}</td></tr>' for a in AF_ITEMS)

# ---- Availability 段全文（去标签，用于报告展示）----
i_av = ms.find("<h3>Availability of data and materials</h3>")
i_av_end = ms.find("<h3>Competing interests</h3>")
AV_TEXT = ms[i_av:i_av_end] if i_av >= 0 else ""
AV_TEXT = re.sub(r"<br\s*/?>", "\n", AV_TEXT)
AV_TEXT = re.sub(r"<[^>]+>", "", AV_TEXT)
AV_TEXT = AV_TEXT.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
AV_TEXT = re.sub(r"\n{3,}", "\n\n", AV_TEXT).strip()
AV_TEXT = AV_TEXT.replace("Availability of data and materials\n", "", 1)
AV_N_GEO = len(re.findall(r"acc\.cgi\?acc=GSE\d+", ms[i_av:i_av_end]))

CSS = """
:root{--blue:#2563eb;--teal:#0d9488;--amber:#d97706;--rose:#e11d48;
--violet:#7c3aed;--green:#16a34a;--ink:#1e293b;--mut:#64748b;--bg:#f8fafc;--line:#e2e8f0}
*{box-sizing:border-box}
body{font-family:"Segoe UI","Microsoft YaHei",sans-serif;background:var(--bg);color:var(--ink);
margin:0;padding:24px;line-height:1.75;font-size:15px}
.wrap{max-width:1000px;margin:0 auto}
h1{font-size:23px;border-bottom:3px solid var(--blue);padding-bottom:10px;margin:0 0 6px}
h2{font-size:18px;color:var(--blue);margin:34px 0 10px}
h3{font-size:15.5px;color:var(--teal);margin:22px 0 8px}
.small{font-size:12.5px;color:var(--mut)}
.box{border:1px solid var(--line);border-left:5px solid var(--blue);background:#fff;border-radius:8px;
padding:14px 18px;margin:14px 0}
.box.def{border-left-color:var(--blue);background:#eff6ff}
.box.tip{border-left-color:var(--teal)}
.box.warn{border-left-color:var(--rose);background:#fff5f5}
.box.rule{border-left-color:var(--amber);background:#fffbeb}
.box.plan{border-left-color:var(--violet);background:#f5f3ff}
.box.ok{border-left-color:var(--green);background:#f0fdf4}
.tag{display:inline-block;background:var(--blue);color:#fff;font-size:11.5px;padding:2px 10px;
border-radius:10px;margin-right:8px;vertical-align:2px}
.box.tip .tag{background:var(--teal)} .box.warn .tag{background:var(--rose)}
.box.rule .tag{background:var(--amber)} .box.plan .tag{background:var(--violet)}
.box.ok .tag{background:var(--green)} .box.def .tag{background:var(--blue)}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.2px}
th{background:#1e293b;color:#fff;padding:8px 10px;text-align:left;font-weight:600}
td{border:1px solid var(--line);padding:7px 10px;background:#fff;vertical-align:top}
tr:nth-child(even) td{background:#f8fafc}
td.c{text-align:center}
.cit{font-size:12.4px;line-height:1.5}
.mono{font-family:Consolas,monospace;font-size:12.8px;background:#f1f5f9;padding:1px 6px;border-radius:4px}
pre{white-space:pre-wrap;word-break:break-word;overflow-wrap:anywhere;background:#fff;
border:1px solid var(--line);border-radius:6px;padding:10px 12px;margin:8px 0;
font-family:Consolas,"Courier New",monospace;font-size:12.2px;line-height:1.62}
.chk{font-size:0.85em;font-weight:700;border-radius:4px;padding:1px 7px;color:#fff}
.chk.y{background:var(--green)} .chk.n{background:var(--rose)} .chk.w{background:var(--amber)}
blockquote{margin:8px 0;padding:6px 14px;border-left:3px solid var(--line);color:#334155;
font-size:13.2px;background:#fff}
.foot{margin-top:36px;border-top:1px solid var(--line);padding-top:12px;font-size:12.5px;color:var(--mut)}
@media print{
  body{background:#fff;padding:0;font-size:11.5px}
  .box,table,tr,blockquote{break-inside:avoid}
  h1,h2,h3{break-after:avoid}
  thead{display:table-header-group}
  @page{margin:14mm}
}
"""

BODY = f"""
<h1>BMC Genomics 投稿准备报告</h1>
<p class="small">稿件：<i>Composition or cell state? Decomposing the adenoma&ndash;carcinoma
transcriptome across tissue, epithelium and single cells</i>｜目标期刊 <b>BMC Genomics</b>（Springer Nature，SCIE / 中科院 2 区 TOP）<br>
本报告依据 BMC Genomics 官网 Submission guidelines 逐条核对｜2026-10-16｜小波</p>

<div class="box ok"><span class="tag">结论</span>
<b>稿件已按 BMC Genomics 要求改造完成，可进入投稿准备。</b>
本轮把一份"内容已达发表标准、但格式仍是内部工作稿"的文件，转成了符合 BioMed Central
规范的 Research article：摘要改三段式并压到 350 词以内、正文引用从哈佛 author-year
整体迁移为 Vancouver 数字制（21 条参考文献按首次出现顺序重新编号）、补齐 BMC 必备的
Conclusions / List of abbreviations / Declarations / Additional files、清掉所有撰写脚手架标签、
图号按正文首见序重排。<br>
<b>随后又补完两处此前留空的部分</b>：① Data availability 从"向通讯作者索取"改写为
完整的<b>可复现链</b>（10 个 GEO 永久标识 + 直链 + 访问日期 + 按 BMC 软件字段格式列出的代码块）；
② Additional files 命名由"正文一套、上传一套"统一为 <span class="mono">Additional file 1&ndash;8</span>，
并按<b>正文首见序</b>自动编号、补齐四要素。自检脚本同步扩到 <b>21 项</b>，全绿。<br>
<b>代码仓库已上线</b>：<span class="mono">github.com/Yao204-116/crc-adenoma-carcinoma-compartment-attribution</span>
（公开 · MIT · 91 文件），Code availability 中的仓库地址已回填为真实值。<br>
<b>剩余待办只有一类：作者信息类字段（姓名、单位、基金、CRediT、ORCID）与 Zenodo 归档 DOI 需本人填写。</b>
</div>

<h2>一、BMC 官方硬指标逐条对照</h2>
<p class="small">下表右列为 BMC Genomics 官网原文核对结果，均引自 <span class="mono">bmcgenomics.biomedcentral.com/submission-guidelines</span>。</p>
<table>
<tr><th style="width:22%">项目</th><th style="width:41%">BMC Genomics 要求（官网原文）</th><th style="width:37%">本稿状态</th></tr>
<tr><td>摘要长度与结构</td>
<td>"The Abstract should not exceed 350 words… The abstract must include the following separate sections: <b>Background</b>: the context and purpose of the study; <b>Results</b>: the main findings; <b>Conclusions</b>: a brief summary and potential implications"</td>
<td><span class="chk y">合规</span> 三段式（无 Methods 小节），<b>{ABS_WORDS} 词</b>，余量 {350-ABS_WORDS} 词</td></tr>
<tr><td>摘要内引用</td><td>"do not cite references in the abstract"</td>
<td><span class="chk y">合规</span> 摘要内零文献角标（原有 GSE 号也已移出）</td></tr>
<tr><td>关键词</td><td>"Three to ten keywords"</td>
<td><span class="chk y">合规</span> 8 个（原 11 个已删减）</td></tr>
<tr><td>主文节序</td><td>Title page, Abstract, Keywords, <b>Background, Methods, Results, Discussion, Conclusions</b>, List of abbreviations, Declarations, References, Figure legends, Tables, Figures</td>
<td><span class="chk y">合规</span> 11 节顺序全对；新增独立 Conclusions（原 Discussion 末段拆出）</td></tr>
<tr><td>参考文献制式</td><td>"References should be cited in the text using <b>consecutive numbers in square brackets</b>. Please use <b>Vancouver</b> reference style."</td>
<td><span class="chk y">合规</span> 21 条，方括号编号，按首次出现顺序；原为哈佛 author-year</td></tr>
<tr><td>作者列名规则</td><td>"If there are more than six authors, list the first six authors followed by et al."</td>
<td><span class="chk y">合规</span> 已用 E-utilities 重取全部 21 篇作者名单，>6 位者列前 6 位 + et al.</td></tr>
<tr><td>图号顺序</td><td>"Figures should be numbered in the order they are first mentioned in the text, and uploaded in this order."</td>
<td><span class="chk y">合规</span> 原顺序 1→6→3→4→5→2→7→8；已重排为严格递增</td></tr>
<tr><td>图题 / 图注字数</td><td>"Figure titles (max 15 words) and legends (max 300 words)"</td>
<td><span class="chk y">合规</span> 11 条图注 57&ndash;161 词；图题最长 9 词</td></tr>
<tr><td>图片文件格式</td><td>接受 EPS / PDF / Word / PPT / TIFF / JPEG / <b>PNG</b> / BMP / CDX；"image resolution of approximately 300 dpi at the final size"；单文件 &le; 10 MB</td>
<td><span class="chk w">待上传</span> 现有 11 个 PNG 均为 300 dpi。PNG 被官方接受，<b>不必强行转 TIFF</b></td></tr>
<tr><td>Declarations 小节</td><td>必备 8 小节：Ethics approval and consent to participate / Consent for publication / Availability of data and materials / Competing interests / Funding / Authors' contributions / Acknowledgements / Authors' information</td>
<td><span class="chk y">合规</span> 8 小节齐全，文字已按"纯公共数据再分析"场景写好</td></tr>
<tr><td>主文词数</td><td>Research article 无词数上限</td>
<td><span class="chk y">合规</span> 主文 {MAIN_WORDS:,} 词，<b>无需瘦身</b></td></tr>
<tr><td>数据共享</td><td>"submission to a community-endorsed, public repository is mandatory"（基因组领域）；数据集须带持久标识与 https 链接；代码须有仓库地址与归档 DOI</td>
<td><span class="chk y">合规</span> 本稿无新数据；10 个公共数据集已给 accession + 直链；衍生表作 Additional files；代码块已写好，仅剩仓库 URL / Zenodo DOI / 许可 3 处占位（见第五节）</td></tr>
</table>

<h2>二、本轮实际改动清单（11 项）</h2>
<table>
<tr><th style="width:4%">#</th><th style="width:26%">改动</th><th>具体内容</th></tr>
<tr><td class="c">1</td><td>摘要重写</td>
<td>原四段式（Background / Methods / Results / Conclusions）431 词压缩为 BMC 三段式 <b>{ABS_WORDS} 词</b>；删去摘要内的 GSE 编号与全部缩写符号（P/U/T 已在正文 List of abbreviations 中定义）；四个关键数字（31.6% / 2.04× / 6&frasl;6 / AUC 0.111）全部保留</td></tr>
<tr><td class="c">2</td><td>引用制式迁移</td>
<td>哈佛 author-year &rarr; Vancouver 数字制：<b>16 处</b>圆括号引用改为 <span class="mono">[n]</span>（连续编号自动压成区间，如 <span class="mono">[1&ndash;4]</span>）；<b>13 处</b>叙事式引用（"X and colleagues"）改为 "X et al. [n]"；Methods 中 WGCNA / Enrichr / MSigDB 三处补软件引用角标（[17][18, 19][20]）</td></tr>
<tr><td class="c">3</td><td>参考文献重排</td>
<td>由字母序 17 篇 + 软件 4 篇改为<b>按正文首次出现顺序编号的 21 条</b>；全部经 E-utilities 重取完整作者名单（原先多篇截断为 3 位作者 + et al.，BMC 规则是 6 位）；期刊缩写改用 PubMed 官方缩写（如 <span class="mono">J Pathol</span> / <span class="mono">Nat Commun</span>）；按官方示例去掉期刊文章 DOI</td></tr>
<tr><td class="c">4</td><td>补独立 Conclusions</td>
<td>原 Discussion 最后一段（"Paragraph 6 &mdash; Conclusion"）拆为独立 <b>Conclusions</b> 章节 &mdash; BMC 要求 Discussion 与 Conclusions 分开</td></tr>
<tr><td class="c">5</td><td>补 List of abbreviations</td>
<td>新增 24 条缩写表（AUC / BH / CAF / DEFA5 / EMT / FAP / LCM / WGCNA / P / U / T 等），BMC 生物学期刊强制</td></tr>
<tr><td class="c">6</td><td>补 Declarations</td>
<td>新增 8 个必备小节（详见上表 3.1）；伦理一栏按"纯公共去标识数据二次分析、无需额外审批"写明，不留空</td></tr>
<tr><td class="c">7</td><td>清撰写脚手架</td>
<td>清掉正文中残留的 <b>9 处</b> <span class="mono">&lt;b&gt;Paragraph N &mdash; …&lt;/b&gt;</span>、<b>11 处</b> <span class="mono">MS ¶N</span> 标签、<b>6 处</b> <span class="mono">¶N</span> 标记，以及 <b>20 处</b>内部章节代码（<span class="mono">（R1）&ndash;（R5）</span>、<span class="mono">M1&ndash;M12</span>、<span class="mono">see Methods M10</span> 等）&mdash; 这些是撰写期的分节标记，投稿稿里不能出现</td></tr>
<tr><td class="c">8</td><td>Methods 小节标题去码</td>
<td>12 个 Methods 小节标题去掉 <span class="mono">M1. &hellip;</span> 前缀，改为短小标题（BMC 要求 "short, informative headings"）</td></tr>
<tr><td class="c">9</td><td>图号重排</td>
<td>Fig. 2 与 Fig. 6 互换（详见第四节），正文引用与图注顺序同步更新</td></tr>
<tr><td class="c">10</td><td>补 Data availability 可复现链</td>
<td>原一句"available from the corresponding author on reasonable request"（BMC <b>明确不接受</b>）扩为：10 个 GEO 数据集逐个带 https 直链 + 访问日期；数据以附加文件形式提供的法定表述；按 BMC 软件字段格式列出的代码块（仓库 / 归档 DOI / OS / 语言 / 依赖 / 许可）。Methods 内一份重复的旧声明同步收口（详见第五节）</td></tr>
<tr><td class="c">11</td><td>Additional files 命名统一 + 按首见序编号</td>
<td>正文 <span class="mono">Fig. S1&ndash;S3 / Supplementary Table S1&ndash;S5</span> 全部改为 <span class="mono">Additional file 1&ndash;8</span>；编号由脚本按<b>正文首见序自动计算</b>（原为"先图后表"，不合规）；清单补齐 BMC 四要素（File name / format / title / description）；附图图注标题改为 <span class="mono">Additional file n (Figure Sn)</span>；并给原 Table S3 补了一处正文引用（详见第六节）</td></tr>
</table>

<h2>三、参考文献编号映射表（21 条）</h2>
<p class="small">编号 = 正文首次出现顺序。全部经 PubMed E-utilities efetch 核实，无编造 PMID。</p>
<table>
<tr><th style="width:5%">#</th><th style="width:58%">BMC / Vancouver 引文</th><th style="width:9%">PMID</th><th style="width:9%">类型</th><th style="width:19%">正文落点</th></tr>
{REF_ROWS}
</table>

<h2>四、图号重排映射</h2>
<div class="box rule"><span class="tag">为什么改</span>
BMC 官网明文："Figures should be numbered in the order they are first mentioned in the text."
按投稿版正文实测的首见顺序为 <b>1 &rarr; 6 &rarr; 3 &rarr; 4 &rarr; 5 &rarr; 2 &rarr; 7 &rarr; 8</b>，
因此只需 <b>Fig. 2 与 Fig. 6 互换</b>，其余主图与全部附图均已合规。改动只涉及 2 处正文引用、
2 条图注标题、2 个图形文件名。
</div>
<table>
<tr><th>原图号</th><th>内容</th><th>正文首见位置</th><th>新图号</th><th>需改动</th></tr>
{FIG_ROWS}
</table>
<div class="box ok"><span class="tag">已同步</span>
<code>figs/</code> 下文件名已按正文图号改到位，并已<b>全部英文化</b>（图内文字与文件名均无中文）：
<span class="mono">fig2_56_pair_concordance.png</span>、
<span class="mono">fig6_three_layer_pattern_adenoma_peak.png</span> 等，与正文图号一致。<br>
另已生成可直接上传的 <span class="mono">upload/</span> 包：主图按 <span class="mono">Figure_1&ndash;8.png</span>
顺序命名（内含 <span class="mono">UPLOAD_MANIFEST.md</span> 上传清单），全部为复制，原文件未动。
</div>

<h2>五、数据可复现链（Availability of data and materials 全文）</h2>
<div class="box rule"><span class="tag">BMC 硬规定</span>
BMC <b>不接受</b>只写"available from the corresponding author on reasonable request"。
数据集须给出 "<i>unique persistent identifier and hyperlink to dataset(s) in https:// format</i>"；
自定义代码须给出仓库地址与<b>带 DOI 的归档版本</b>（官网明确推荐 GitHub 项目用 Zenodo 归档取得 DOI）。
本轮已按此逐条改写，把原来一句"向通讯作者索取"扩成了完整的可复现链。
</div>
<table>
<tr><th style="width:15%">链条环节</th><th style="width:43%">内容</th><th>落在稿件的哪里</th></tr>
<tr><td><b>① 原始数据</b></td>
<td>10 个 GEO 数据集逐个给出 accession + 完整 https 直链（共 <b>{AV_N_GEO}</b> 条）+ 统一访问日期，并按"发现/全组织验证队列、单细胞队列、空间参照"分组</td>
<td>Declarations &rarr; Availability of data and materials（自检已核：10/10 直链齐全）</td></tr>
<tr><td><b>② 衍生数据</b></td>
<td>8 个 additional file（5 张表 + 3 张图），含 BMC 要求的四要素，全部机器可读（.xlsx / .csv）</td>
<td>稿件第 11 节 Additional files（详见下一节映射表）</td></tr>
<tr><td><b>③ 代码</b></td>
<td>仓库地址 + Zenodo 归档 DOI + 操作系统 + 编程语言 + 依赖 + 许可，按 BMC 的软件字段格式逐项列出</td>
<td>同段 <i>Code availability</i> 子块（<span class="chk n">3 处待填</span>）</td></tr>
<tr><td><b>④ 复现清单</b></td>
<td>脚本&rarr;产物的逐条对应、运行环境、外部前置依赖、Zenodo 提交清单、以及"哪些结论不靠代码"的说明</td>
<td><span class="mono">REPRODUCIBILITY.md</span>（本工作区新建，随代码一起入库）</td></tr>
</table>
<div class="box def"><span class="tag">稿件中现在的原文</span>
<pre>{AV_TEXT}</pre>
</div>
<div class="box warn"><span class="tag">作者还需填 3 处</span>
① 仓库地址 <span class="mono">[Author to insert: public repository URL]</span>；
② Zenodo 归档 DOI <span class="mono">[Author to insert: DOI of the frozen release]</span>；
③ 许可 <span class="mono">[Author to choose, e.g. MIT]</span>（R 版本也建议一并确认）。<br>
另有一处<b>容易踩坑</b>：<span class="mono">build_suppl_tables.py</span> 依赖两个本机 R 流程的产物
<span class="mono">D:/2/02.WGCNA/WGCNA_top_module_genes.csv</span> 与
<span class="mono">D:/2/03.enrich/Hub_Tumor_vs_Polyp.csv</span>，
<b>入库时必须一并带上</b>（或把上游 R 脚本一起放进去），否则 Additional file 2 无法重建 ——
审稿人一旦试着复现就会卡在这里。
</div>

<h2>六、Additional files 编号与上传映射</h2>
<div class="box rule"><span class="tag">BMC 两条硬规定</span>
① "Each additional file should be cited <b>in sequence</b> within the main body of text."
&mdash;&mdash; 编号必须等于<b>正文引用顺序</b>；<br>
② "Additional files should be named <b>'Additional file 1'</b> and so on and should be
<b>referenced explicitly by file name</b> within the body"；清单须含四项：
<i>File name / File format / Title of data / Description of data</i>。
</div>
<div class="box warn"><span class="tag">本轮改了什么</span>
<b>改前</b>：正文写 <span class="mono">Fig. S1&ndash;S3</span> /
<span class="mono">Supplementary Table S1&ndash;S5</span>，上传却叫
<span class="mono">Additional file 1&ndash;8</span> &mdash;&mdash; 两套命名并存，审稿人无法对应；
且编号顺序是"先图后表"，与正文首见序不符（正文最早引用的其实是 Table S1）。<br>
<b>改后</b>：正文全部统一为 <span class="mono">Additional file n</span>；编号由脚本按正文首见序
<b>自动计算</b>（避免手工排错）；清单补齐四要素；附图图注标题改写为
<span class="mono">Additional file n (Figure Sn)</span>。<br>
另补了一处正文引用：原 Table S3（三层模式矩阵）此前<b>只在 Fig. 6 图注里出现</b>，
不满足"cited within the main body"，现已在 Results 的三层模式句内引用。
</div>
<table>
<tr><th style="width:9%">编号</th><th style="width:8%">内部标签</th><th style="width:26%">数据标题</th>
<th style="width:8%">格式</th><th style="width:9%">正文落点</th>
<th style="width:18%">上传文件名建议</th><th style="width:22%">工作区源文件</th></tr>
{AF_ROWS}
</table>
<div class="box tip"><span class="tag">注意：AF 号 ≠ S 号</span>
附加文件编号跟随<b>正文引用顺序</b>，因此 Table S4 / S5 排到了 Table S3 之前
（<span class="mono">Additional file 3 = Table S4</span>、
<span class="mono">Additional file 7 = Table S3</span>）。
这是 BMC 规则决定的，不是笔误 &mdash;&mdash; 清单里同时保留内部标签，读者与编辑都能对上。
如果你更希望"AF 号与 S 号一一对应"，那就必须给 Table S3 补一处<b>更靠前</b>的正文引用，
编号才会跟着前移；否则保持现状更稳。
</div>

<h2>七、必须由作者补齐的信息</h2>
<div class="box warn"><span class="tag">投稿系统会卡住</span>
以下字段本稿以方括号占位（稿件首页与 Declarations 内），<b>投稿前必须由本人填写</b>，否则系统无法提交或会被编辑部退回补件：
</div>
<table>
<tr><th style="width:24%">字段</th><th>当前占位</th><th style="width:34%">怎么填</th></tr>
<tr><td>作者姓名与顺序</td><td><span class="mono">[Author 1] [Author 2] [Author 3]…</span></td><td>按最终署名顺序；BMC 会要求每位作者的邮箱</td></tr>
<tr><td>单位与通讯地址</td><td><span class="mono">[Department, Institution, City, Country]</span></td><td>全称 + 邮编；单位内部不同科室用上标区分</td></tr>
<tr><td>通讯作者信息</td><td><span class="mono">*Corresponding author: [name], [e-mail], [ORCID iD]</span></td><td>BMC 投稿系统要求通讯作者 ORCID</td></tr>
<tr><td>基金</td><td>已写"未获任何专项资助"，另有 <span class="mono">[Author to confirm…]</span></td><td>若有基金必须补编号；确实没有则删掉方括号句</td></tr>
<tr><td>作者贡献（CRediT）</td><td><span class="mono">[Initials] Conceptualization; …</span></td><td>用姓名缩写填 7 类角色，导师通常填 Supervision</td></tr>
<tr><td>数据与代码位置</td><td>原为"available from the corresponding author on reasonable request"</td><td><b>已改写</b>为完整可复现链（10 个 GEO 直链 + 代码块），仅剩 3 处占位：仓库地址、Zenodo DOI、许可（详见第五节）</td></tr>
</table>

<h2>八、最高策略风险：BMC 的"纯再分析"排除条款</h2>
<div class="box warn"><span class="tag">必读</span>
BMC Genomics 的 Research article 收稿标准里有一句：<br>
<blockquote>"Please note that non-commissioned <b>pooled analyses of selected published research</b>
and <b>bibliometric analyses</b> will not be considered."</blockquote>
本稿是公共数据再分析，<b>从措辞上有被套进"pooled analysis"的可能</b>。必须提前把界限讲清楚——
否则有 desk reject 风险。<br><br>
<b>我们的实际情况与"pooled analysis"的区别（三点，均可在 Cover Letter 里直说）：</b><br>
① <b>分析对象不同</b>：pooled analysis 汇总的是<b>已发表论文的结果</b>；本稿直接分析
<b>GEO 原始测序/芯片数据</b>（CEL / 表达矩阵），自行做探针折叠、批次处理、漂移校正，不是把别人的表格拼起来。<br>
② <b>有独立的方法学约束</b>：同一签名在三种分辨率下打分，并要求在两个口径（未校正/漂移校正）
下同时复现，另设两个必须"什么都测不出来"的阴性对照——这是<b>事先设定的设计纪律</b>，不是事后筛选文献。<br>
③ <b>有可裁决的具体结论</b>：把腺瘤期分化峰归属到 per-cell state 而非细胞丰度，并指出髓系
<span class="mono">LYZ</span> 会把全细胞 Paneth 读数在腺瘤&rarr;癌之间读反——这是<b>新结论</b>，不是文献综述。<br><br>
<b>另一条保命线</b>：BMC Genomics 官网明文它不按"重要性/趣味性"做编辑决定，只看
"scientifically valid"（研究问题清晰 + 方法学可靠 + 符合领域共识标准）。本稿的优势恰好在这三点上。
</div>

<h2>九、投稿前必办清单</h2>
<table>
<tr><th style="width:4%">#</th><th style="width:30%">事项</th><th>说明</th></tr>
<tr><td class="c">1</td><td><span class="chk n">必办</span> 填齐第七节全部占位字段</td><td>姓名、单位、通讯作者、基金、CRediT</td></tr>
<tr><td class="c">2</td><td><span class="chk y">已办</span> 图文件英文化 + 图号 + 上传包</td><td><b>10 张图的图内文字与文件名已全部英文化</b>（原为中文，投英文刊属阻塞项），图内图号与正文一致；<span class="mono">upload/</span> 已生成 <span class="mono">Figure_1&ndash;8.png</span>（按正文图号排序）+ <span class="mono">Additional_file_1&ndash;8.*</span>。PNG 已被 BMC 接受，若想更稳可另存 EPS/PDF 矢量版</td></tr>
<tr><td class="c">3</td><td><span class="chk y">已办</span> 代码入库（含 DOI 待取）</td><td>仓库已上线：<span class="mono">github.com/Yao204-116/crc-adenoma-carcinoma-compartment-attribution</span>（公开 · MIT · 91 文件），稿件 Code availability 的仓库地址已回填。★ <b>剩一步</b>：在 Zenodo 打开该仓库开关后发 release 取归档 DOI（必须先开开关再发 release），再把 DOI 回填</td></tr>
<tr><td class="c">4</td><td><span class="chk y">已查</span> 投稿前两处文献状态（2026-10-16 实测）</td><td>① <b>GSE303070 关联论文仍是预印本</b>：Pattison AD, …, Abud HE. <i>A spatial atlas of colorectal cancer reveals the influence of stromal niches on tumour differentiation.</i> bioRxiv 2025；DOI <span class="mono">10.64898/2025.12.09.693138</span>（Europe PMC <span class="mono">PPR1231869</span>）——<b>未进入 PubMed，GEO 记录亦无 PMID</b>。文中"nine spatial niches"一说出自该数据集的 GEO summary（原文 "identified 9 spatial niches"），<b>以 GEO accession 归属站得住</b>；若要更规范，可在该句后补一条参考文献（Vancouver 顺序下会成为 <span class="mono">[22]</span>）。<br>② <b>Guha et al. 仍是预印本</b>：PubMed <span class="mono">PMID 39605357</span>，bioRxiv，article type 含 <i>Preprint</i>，Europe PMC 无期刊版 → 现用的 bioRxiv 引用格式（<span class="mono">[8]</span>）<b>正确，无需改动</b>。<br>③ 附带一条：Guha TK 另有 <b>Nature 2026</b> 论文 <i>Polyclonal origins of human premalignant colorectal lesions</i>（<span class="mono">PMID 41291291</span>）——同组不同题，<b>不替代</b>上述预印本，勿误改引用</td></tr>
<tr><td class="c">5</td><td><span class="chk y">已办</span> Additional files 命名统一</td><td>正文已全改为 <span class="mono">Additional file 1&ndash;8</span> 并按首见序编号（映射见第六节）；上传时用 <span class="mono">Additional_file_n.*</span> 命名即可</td></tr>
<tr><td class="c">6</td><td><span class="chk y">已办</span> 摘要 &le;350 词</td><td>{ABS_WORDS} 词；若编辑另有要求，优先压缩 Background 第二句</td></tr>
<tr><td class="c">7</td><td><span class="chk y">已办</span> 引用编号闭环</td><td>21 条全部被正文引用，正文角标全部有对应条目（脚本 <span class="mono">verify_bmc.py</span> 自动核对）</td></tr>
<tr><td class="c">8</td><td><span class="chk y">已办</span> 零孤儿图</td><td>11 张图全部在正文被引用，且首见序严格递增</td></tr>
</table>

<h2>十、Cover Letter 初稿（BMC Genomics 定制）</h2>
<div class="box plan"><span class="tag">可直接改用</span>
方括号处需替换。第二段刻意回应了"pooled analysis"排除条款。
</div>
<blockquote>
Dear Editor,<br><br>
We submit for your consideration our manuscript entitled
"Composition or cell state? Decomposing the adenoma&ndash;carcinoma transcriptome across
tissue, epithelium and single cells"
for publication as a Research article in <i>BMC Genomics</i>.<br><br>
Colorectal carcinogenesis proceeds through the adenoma&ndash;carcinoma sequence, yet expression
differences between these stages are almost always measured in whole tissue, where epithelial,
stromal and immune contributions are pooled. A tissue-level difference can therefore arise from a
change in cellular composition, from a change in the programme executed by the epithelial cells
themselves (state), or from both. This distinction must be resolved before specimens of differing
epithelial content are compared, and at the adenoma stage it has not been. We derive a 313-gene
signature of the transition from four independent cohorts and score it at three resolutions, each
fixing one term of the composition&ndash;state relationship: whole tissue; laser-capture-microdissected
epithelium, used as a physical control in which the stromal term is removed rather than estimated;
and patient-matched single-cell RNA sequencing of normal&ndash;adenoma&ndash;carcinoma triads. The two arms
of the signature have opposite cellular sources, and the adenoma-stage differentiation peak is
epithelial-intrinsic and carried by per-cell state rather than by cell abundance. We also audited an
independent spatial atlas and report it as a negative result: such a resource can mislead as a
composition benchmark when its sampling is not balanced across stages.<br><br>
We wish to make the nature of the analysis explicit, since the journal's criteria exclude
non-commissioned pooled analyses of previously published research. This study is not such an
analysis. It does not aggregate or re-analyse published results; it processes the underlying
primary measurements deposited in the NCBI Gene Expression Omnibus, applying its own probe
collapse, cross-platform alignment and drift correction. Its design constraints were fixed in
advance rather than selected after the fact: every panel trajectory had to reproduce in all three
cohorts under two independent scoring calibers, and two panels designed to be uninformative served
as negative controls. The work reports a specific, falsifiable claim rather than a literature
synthesis &mdash; the quantitative attribution of an established progression signal to composition
versus state at a single transition, together with the identification of myeloid <i>LYZ</i> as an
artefact that inverts the whole-cell Paneth readout between adenoma and carcinoma.<br><br>
The study fits <i>BMC Genomics</i> because it is a genome-scale transcriptomic analysis whose
conclusions rest on cross-platform reproducibility and explicit negative controls, and because the
journal evaluates manuscripts on scientific validity rather than perceived impact. Every analysis
uses publicly available, de-identified data; no new data were generated. The manuscript has not been
published previously and is not under consideration elsewhere. The authors declare no competing
interests.<br><br>
Thank you for your consideration.<br><br>
Sincerely,<br>
[Corresponding author], on behalf of all authors<br>
[Affiliation &middot; e-mail &middot; ORCID]
</blockquote>

<h2>十一、可复现脚本</h2>
<table>
<tr><th style="width:34%">脚本</th><th>作用</th></tr>
<tr><td><span class="mono">_fetch_bmc_refs.py</span></td><td>E-utilities efetch 重取 21 篇完整作者名单，按 BMC Vancouver 格式生成 <span class="mono">_bmc_references.json</span></td></tr>
<tr><td><span class="mono">build_bmc_manuscript.py</span></td><td>把《全文整合_投稿版》改造为 <span class="mono">Manuscript_BMC_Genomics_投稿版.html</span>（引用迁移 / 摘要重写 / 节序 / 脚手架清理 / 图号重排 / <b>附加文件按首见序自动编号</b> / <b>Availability 可复现链</b>）</td></tr>
<tr><td><span class="mono">verify_bmc.py</span></td><td><b>21 项</b>合规自检：残留 author-year、汉字、脚手架、编号闭环、<b>引文按正文首见序编号（Vancouver）</b>、图号顺序、摘要与图注字数、节序、Declarations 完整性、<b>附加文件首见序与四要素、Availability 的 10 个 GEO 直链完备性</b></td></tr>
<tr><td><span class="mono">REPRODUCIBILITY.md</span></td><td>随代码入库的复现清单：环境 / 10 个公共数据集 / 9 个分析阶段（脚本&rarr;产物逐条对应）/ 单条产物复现示例 / 外部前置依赖 / Zenodo 提交清单</td></tr>
</table>

<div class="foot">
生成于 2026-10-16｜依据 BMC Genomics 官网 Submission guidelines（2026-10 核对）｜小波<br>
配套文件：<span class="mono">Manuscript_BMC_Genomics_投稿版.html</span>（投稿稿）&middot;
<span class="mono">Manuscript_BMC_Genomics_投稿版.pdf</span>（桌面）&middot;
<span class="mono">BMC投稿准备报告_BMC_Genomics.pdf</span>（桌面）
</div>
"""

out = f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<title>BMC Genomics 投稿准备报告</title>
<style>{CSS}</style></head><body><div class="wrap">
{BODY}
</div></body></html>"""

open(OUT, "w", encoding="utf-8").write(out)
print(f"wrote {OUT} ({os.path.getsize(OUT):,} bytes)")
print(f"abstract words={ABS_WORDS}  main words={MAIN_WORDS}  refs={len(refs)}")
for r in refs:
    print(f"  {r['num']:2d}. {where(r['num'])}")
