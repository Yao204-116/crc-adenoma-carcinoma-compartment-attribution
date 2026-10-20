# -*- coding: utf-8 -*-
"""
把《全文整合_投稿版》改造成 BMC Genomics Research article 投稿稿。

依据（BMC Genomics 官方 Submission guidelines，逐字）：
  * Abstract: "should not exceed 350 words" + "must include the following separate
    sections: Background / Results / Conclusions" + "do not cite references in the abstract"
  * Keywords: "Three to ten keywords"
  * 主文节序: Title page, Abstract, Keywords, Background, Methods, Results, Discussion,
    Conclusions, List of abbreviations, Declarations, References, Figure legends, Tables
  * References: "cited in the text using consecutive numbers in square brackets.
    Please use Vancouver reference style."
  * Figures: "should be numbered in the order they are first mentioned in the text,
    and uploaded in this order"; "Figure titles (max 15 words) and legends (max 300 words)"
  * Declarations 必备 8 小节

本脚本做 9 件事：
  1. 清掉撰写脚手架标签（<b>Paragraph N — …</b> / <span class="tag">MS ¶N</span> / <span class="lbl">¶N</span>）
  2. 摘要换成 BMC 三段式（≤350 词、无文献引用）
  3. Introduction → Background
  4. 正文引用：哈佛 author-year → Vancouver 方括号编号
  5. 补软件引用（WGCNA / Enrichr / MSigDB）与叙事式引用角标
  6. Results 补 3 条小节标题；Discussion ¶6 拆出为独立 Conclusions
  7. 新增 List of abbreviations / Declarations / Additional files
  8. References 换成按引用顺序编号的 Vancouver 列表（21 条）
  9. 图号按正文首见序重排（Fig. 2 ↔ Fig. 6）
"""
import re, json, os

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)

SRC = "Manuscript_全文整合_投稿版.html"
OUT = "Manuscript_BMC_Genomics_投稿版.html"

# ---------------------------------------------------------------- 引用编号表
# 编号 = 正文首次出现顺序（Background → Methods → Results）
CITE_NUM = {
    ("McCorry", "2018"): 1, ("Kreis", "2024"): 2, ("Tyler", "2021"): 3,
    ("Block", "2019"): 4, ("Yan", "2022"): 5, ("Becker", "2022"): 6,
    ("Zheng", "2022"): 7, ("Guha", "2025"): 8, ("Yang", "2025"): 9,
    ("Bala", "2023"): 10, ("Porter", "1997"): 11, ("Mallow", "1996"): 12,
    ("Radeva", "2010"): 13, ("Joo", "2009"): 14, ("Pagnini", "2011"): 15,
    ("López-Arribillaga", "2021"): 16, ("Langfelder", "2008"): 17,
    ("Kuleshov", "2016"): 18, ("Chen", "2013"): 19, ("Liberzon", "2015"): 20,
    ("Thiruvengadam", "2019"): 21,
}

# 叙事式引用（"X and colleagues" / "Tyler and Tirosh"）→ 编号
NARRATIVE = {
    "McCorry": 1, "Kreis": 2, "Block": 4, "Becker": 6, "Radeva": 13,
    "Joo": 14, "Pagnini": 15, "López-Arribillaga": 16,
    "Thiruvengadam": 21, "Yan": 5, "Yang": 9,
    # 2026-10-20：Discussion ¶3/¶5 引用 GSE109308 时加入
    "Bala": 10,
}

# ---------------------------------------------------------------- 摘要（BMC 三段式）
ABSTRACT = {
    "Background": (
        "Colorectal carcinogenesis proceeds through the adenoma-carcinoma sequence, but "
        "expression differences between these stages are usually measured in whole tissue, "
        "where epithelial, stromal and immune contributions are pooled. Whether a "
        "stage-associated programme reflects altered cellular composition or altered "
        "epithelial cell state has remained unresolved for the adenoma stage."
    ),
    "Results": (
        "We derived a 313-gene signature of the adenoma-to-carcinoma transition from four "
        "independent cohorts and scored it at three resolutions, each fixing one term of the "
        "composition-state relationship: "
        "whole tissue; laser-capture-microdissected epithelium (6 adenomas, 6 carcinomas, "
        "6 normals); and patient-matched single-cell RNA sequencing of three "
        "normal-adenoma-carcinoma triads (13 libraries, 54,782 cells). Panel trajectories "
        "had to reproduce in all three cohorts under both uncorrected and drift-corrected "
        "scoring; two designs that had to yield nothing served as negative controls. The "
        "two components of the signature behaved oppositely: in purified epithelium the "
        "up-regulated component retained only 31.6% of its whole-tissue amplitude (gene "
        "concordance 58.6%, near chance), whereas the down-regulated component was "
        "amplified 2.04-fold (93.2% concordant), assigning the tissue-level increase to "
        "stromal and myeloid composition and the decrease to an epithelial-intrinsic loss "
        "of differentiation. Of thirteen cell-programme panels, only the Paneth-cell "
        "programme (DEFA5, DEFA6, OLFM4, SPINK1) peaked in adenomas in all three cohorts "
        "under both conventions (normal &lt; carcinoma &lt; adenoma) and persisted in "
        "epithelium-only samples. Single-cell decomposition attributed that peak mainly to "
        "per-cell state (58% of the summed component change in one patient; the composition "
        "term was negative in the other), with DEFA6 higher in adenoma than in both the "
        "matched normal and the matched carcinoma in 6/6 paired comparisons. In two patients "
        "the whole-cell Paneth score rose from adenoma to carcinoma while the "
        "intra-epithelial score did not, an artefact of myeloid LYZ entering the whole-cell "
        "average. The negative-control panel gave a null paired effect (28/56, P = 0.70)."
    ),
    "Conclusions": (
        "A progression signature containing oppositely directed components must be "
        "decomposed before it is compared across specimens differing in epithelial content: "
        "in epithelium-only samples the mixed score reversed (AUC 0.111) rather than lost "
        "signal. Decomposed, the adenoma-stage differentiation peak is epithelial-intrinsic "
        "and carried by per-cell state rather than cell abundance, whereas the tissue-level "
        "progression signal is largely compositional."
    ),
}

KEYWORDS = [
    "colorectal adenoma-carcinoma sequence", "laser-capture microdissection",
    "single-cell RNA sequencing", "cell composition", "cell state",
    "Paneth cell", "defensin", "transcriptional signature",
]

ABBREVIATIONS = [
    ("AUC", "area under the receiver operating characteristic curve"),
    ("BH", "Benjamini-Hochberg"),
    ("CAF", "cancer-associated fibroblast"),
    ("CI", "confidence interval"),
    ("DEFA5/DEFA6", "defensin alpha 5 / defensin alpha 6"),
    ("DIA", "data-independent acquisition"),
    ("EMT", "epithelial-mesenchymal transition"),
    ("FAP", "familial adenomatous polyposis"),
    ("FDR", "false discovery rate"),
    ("FOV", "field of view"),
    ("GEO", "Gene Expression Omnibus"),
    ("GO", "Gene Ontology"),
    ("GSE", "GEO Series accession"),
    ("LCM", "laser-capture microdissection"),
    ("MSigDB", "Molecular Signatures Database"),
    ("OLFM4", "olfactomedin 4"),
    ("P", "composition term (cell-abundance component of a mixed signature score)"),
    ("QC", "quality control"),
    ("scRNA-seq", "single-cell RNA sequencing"),
    ("SPINK1", "serine peptidase inhibitor Kazal type 1"),
    ("SSA", "sessile serrated adenoma"),
    ("T", "whole-cell (total) score"),
    ("U", "cell-intrinsic state term"),
    ("WGCNA", "weighted gene co-expression network analysis"),
]

# ---------------------------------------------------------------- Availability 语句
# BMC 要求（逐字）：
#  · 数据集须给出 "unique persistent identifier and hyperlink to dataset(s) in https:// format"
#  · 数据以附加文件形式提供时，须写 "The dataset(s) supporting the conclusions of this
#    article is(are) included within the article (and its additional file(s))."
#  · 自定义代码须给出 Project home page / Archived version(DOI) / OS / 语言 / 依赖 / 许可
#    （BMC 推荐 GitHub 代码用 Zenodo 归档取得 DOI）
#  · 不允许只写 "available from the corresponding author on reasonable request"
GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
GEO_ACCESSED = "16 October 2026"

# ---- Code availability 的单一数据源 -------------------------------------
# 仓库 URL / Zenodo DOI / 许可 / R 版本只在 `code_availability.json` 里维护一份，
# 由 `push_github.py` 推送成功并取到 DOI 后自动写入；文件不存在时回退为占位文本。
# 这样"重建稿件"永远不会把已填好的 DOI 覆盖回占位符。
CODE_AVAIL_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "code_availability.json")
_PLACEHOLDER = {
    "repo_url": "[Author to insert: public repository URL, "
                "e.g. https://github.com/&lt;user&gt;/&lt;repo&gt;]",
    "doi": "[Author to insert: DOI of the frozen release archived at Zenodo, "
           "e.g. https://doi.org/10.5281/zenodo.XXXXXXX]",
    "license": "[Author to choose, e.g. MIT]",
    "r_version": "[Author to confirm the R version]",
}


def _code_avail():
    vals = dict(_PLACEHOLDER)
    try:
        with open(CODE_AVAIL_FILE, encoding="utf-8") as fh:
            for k, v in json.load(fh).items():
                if v and k in vals:
                    vals[k] = v
    except Exception:
        pass
    return vals


CODE_AVAIL = _code_avail()

_GEO_GROUPS = [
    ("Discovery and whole-tissue validation cohorts",
     ["GSE41258", "GSE37364", "GSE20916", "GSE41657",
      "GSE4183", "GSE117606", "GSE15960", "GSE128435"]),
    ("Single-cell cohorts", ["GSE161277", "GSE109308"]),
    ("Spatial-transcriptomics reference (used only as a descriptive counter-example)",
     ["GSE303070"]),
]


def _availability():
    p = []
    p.append(
        "All datasets analysed in this study are publicly available in the NCBI Gene "
        "Expression Omnibus (GEO) repository and are listed below with their permanent "
        "identifiers. No new sequencing data were generated. The datasets supporting the "
        "conclusions of this article are included within the article and its additional "
        "files (Additional files 1-8), which provide the derived consensus signature gene "
        "list, the per-gene attribution table, the panel-by-cohort discrimination matrix, "
        "the three-layer pattern matrix and the single-cell composition/state decomposition "
        "in machine-readable (.xlsx / .csv) format. Analysis code is available as described "
        "under <i>Code availability</i> below.")
    for title, accs in _GEO_GROUPS:
        items = "; ".join(f"{a} ({GEO_BASE}{a})" for a in accs)
        p.append(f"<b>{title}:</b> {items}.")
    p.append(f"All URLs were last accessed on {GEO_ACCESSED}. No access restrictions or "
             "registration are required for any of these datasets.")
    ca = CODE_AVAIL

    def _link(v):
        """真 URL 给出可点击的 https 链接；占位文本原样保留。"""
        return v if v.lstrip().startswith("[Author") else f'<a href="{v}">{v}</a>'

    p.append("<b>Code availability.</b> "
             "Project name: adenoma-carcinoma compartment attribution pipeline. "
             f"<br>Project home page: {_link(ca['repo_url'])}. "
             f"<br>Archived version: {_link(ca['doi'])}. "
             "<br>Operating system(s): platform independent (developed on Windows 11). "
             "Programming language: Python 3.13 and R (WGCNA and limma packages; "
             f"{ca['r_version']}). "
             "<br>Other requirements: see the environment file in the repository. "
             f"<br>License: {ca['license']}. "
             "<br>Any restrictions to use by non-academics: none.")
    return "\n".join(p)


AVAILABILITY = _availability()


# ---------------------------------------------------------------- 作者元数据（单一数据源）
# ★ 一个文件驱动三处，避免"稿件署名 ≠ DOI 记录署名"：
#     稿件标题页 / Funding / Authors' contributions  ← 本脚本读
#     `.zenodo.json` 的 creators                      ← fill_author_metadata.py 读
#   文件不存在时退回方括号占位（早期通读用，但**绝不能带着它发 release**）。
AUTHOR_META_PATH = os.path.join(BASE, "author_metadata.json")


def _author_meta():
    if not os.path.exists(AUTHOR_META_PATH):
        return None
    try:
        m = json.load(open(AUTHOR_META_PATH, encoding="utf-8"))
    except Exception:
        return None
    if not m.get("authors"):
        return None
    # ★ 模板（姓名还空着）也当"没有"处理 —— 否则标题页会渲染出空名字。
    if not any((a.get("name") or "").strip() for a in m["authors"]):
        return None
    return m


META = _author_meta()

# CRediT 固定排序（BMC 无强制顺序，但列出时按此更易读）
_CREDIT_ORDER = [
    "Conceptualization", "Methodology", "Software", "Validation", "Formal analysis",
    "Investigation", "Resources", "Data curation", "Writing - original draft",
    "Writing - review and editing", "Visualization", "Supervision",
    "Project administration", "Funding acquisition",
]


def _corr_authors(m):
    """全部通讯作者（保持署名顺序）。BMC 允许 co-corresponding，故返回列表不是单值。"""
    cs = [a for a in m["authors"] if a.get("corresponding")]
    return cs or [m["authors"][0]]


def _author_email(m, a, n_corr):
    """作者邮箱：作者自带 email 优先（多通讯写法），否则退回顶层 corresponding.email。"""
    mail = (a.get("email") or "").strip()
    if not mail and n_corr == 1:
        mail = ((m.get("corresponding") or {}).get("email") or "").strip()
    # ★ 缺失时降级为可见占位 `[e-mail]`（而不是渲染成空、也不是 KeyError）——
    #   闸门的"任意非数字方括号"扫描会把它抓出来，所以不会带着它发出去。
    return mail or "[e-mail]"


def title_page_html():
    if not META:
        return (
            '<p class="authline">[Author 1]<sup>1</sup>, [Author 2]<sup>1</sup>, '
            '[Author 3]<sup>2</sup>, [Corresponding author]<sup>1*</sup></p>\n'
            '<p class="authline"><sup>1</sup>[Department, Institution, City, Country] &nbsp; '
            '<sup>2</sup>[Department, Institution, City, Country]</p>\n'
            '<p class="authline">*Corresponding author: [name], [e-mail], '
            '[ORCID iD]</p>\n')
    m = META
    corr = _corr_authors(m)
    n_corr = len(corr)
    names = []
    for a in m["authors"]:
        inner = ",".join(str(n) for n in (a.get("aff") or [1]))
        if a.get("corresponding"):
            inner += "*"          # ★ 通讯符加在 <sup> 的**内容**里，不能切字符串拼标签
        names.append("%s<sup>%s</sup>" % (a["name"], inner))
    affs = ["<sup>%s</sup>%s" % (k, m["affiliations"][k])
            for k in sorted(m["affiliations"], key=lambda x: int(x))]
    # 多通讯：联系人**并列在同一段**（`; ` 分隔），而不是一人一行。
    # 两人都是通讯时，重复 "*Corresponding authors:" 标签会读成语病（"作者们：A。/ 作者们：B。"）。
    corr_bits = []
    for ca in corr:
        s = "%s, %s" % (ca["name"], _author_email(m, ca, n_corr))
        if ca.get("orcid"):
            s += ", ORCID iD: %s" % ca["orcid"]
        corr_bits.append(s)
    if n_corr == 1:
        # 单通讯：沿用旧形态 "*Corresponding author: Name, mail"（首字符即通讯符）。
        # ★ 通讯符**不能**只在 `names` 那行的 `<sup>` 里加 —— 历史上本行是独立渲染的，
        #   单通讯时它自带头一个 `*`；`_corr_authors()` 在"无一人标 corresponding"时会
        #   退回第一位作者，那种情况下 `names` 行也拿不到 `*`，只能靠这里兜底。
        corr_line = "*Corresponding author: " + corr_bits[0]
    else:
        corr_line = "*Corresponding authors: " + "; ".join(corr_bits)
    return "\n".join([
        '<p class="authline">%s</p>' % ", ".join(names),
        '<p class="authline">%s</p>' % " &nbsp; ".join(affs),
        '<p class="authline">%s</p>' % corr_line]) + "\n"


TITLE_META_HTML = (
    '<p class="meta">Research article submitted to <b>BMC Genomics</b> &middot; '
    'Vancouver reference style &middot; 21 references &middot; 8 main figures &middot; '
    '8 additional files (3 supplementary figures, 5 supplementary tables)</p>')


def funding_text():
    if META and META.get("funding"):
        return META["funding"]
    return ("This research did not receive any specific grant from funding agencies in the "
            "public, commercial, or not-for-profit sectors. [Author to confirm: insert grant "
            "numbers or state that no funding was received.]")


def contributions_text():
    """按 CRediT 汇总：{initials: [roles]} → 'Role: A, B; Role2: C'。"""
    if META and META.get("credit"):
        buckets = {}
        for ini, roles in META["credit"].items():
            for r in roles:
                buckets.setdefault(r, []).append(ini)
        seq = [r for r in _CREDIT_ORDER if r in buckets]
        seq += sorted(set(buckets) - set(seq))
        body = "; ".join("%s: %s" % (r, ", ".join(sorted(buckets[r]))) for r in seq)
        return body + ". All authors read and approved the final manuscript."
    return ("[Initials] Conceptualization; [initials] Methodology; [initials] Formal analysis; "
            "[initials] Data curation; [initials] Writing - original draft; [initials] Writing - "
            "review and editing; [initials] Supervision. All authors read and approved the final "
            "manuscript. [Author to complete per CRediT.]")


# BMC Declarations 必备 8 小节；方括号内容由作者填写
DECLARATIONS = [
    ("Ethics approval and consent to participate",
     "Not applicable. This study is a secondary analysis of publicly available, "
     "de-identified human transcriptomic data deposited in the NCBI Gene Expression "
     "Omnibus. No human participants, human material or animals were involved in this "
     "work, and no new data were generated; ethics approval and informed consent were "
     "obtained by the investigators of the original studies."),
    ("Consent for publication",
     "Not applicable. This manuscript contains no individual person's data."),
    ("Availability of data and materials", AVAILABILITY),
    ("Competing interests",
     "The authors declare that they have no competing interests."),
    ("Funding", funding_text()),
    ("Authors' contributions", contributions_text()),
    ("Acknowledgements",
     "The authors gratefully acknowledge the investigators of the public datasets analysed "
     "here (GSE41258, GSE37364, GSE20916, GSE41657, GSE4183, GSE117606, GSE15960, "
     "GSE128435, GSE161277, GSE109308 and GSE303070) for making their data openly "
     "available."),
    ("Authors' information",
     "Not applicable."),
]

# BMC 要求（逐字）："Each additional file should be cited in sequence within the main body
# of text." 且附加文件统一命名 "Additional file 1"、正文须"referenced explicitly by file
# name"。清单须含四项：File name / File format(带扩展名) / Title of data / Description of data。
# 因此编号不写死，改由"正文首见序"自动计算（见 unify_additional_refs）。
#   key, 正文标记, 内部标签, 格式, 标题, 描述
ADDITIONAL_ITEMS = [
    ("T1", "Supplementary Table S1", "Table S1", ".xlsx, .csv",
     "Per-gene attribution to the epithelial and stromal compartments.",
     "One row per consensus gene (313 genes): whole-tissue and LCM epithelial effects, "
     "drift-corrected z-scores, component assignment and direction."),
    ("T2", "Supplementary Table S2", "Table S2", ".xlsx, .csv",
     "Gene-set membership of the panels and the two consensus components.",
     "Gene lists for the 17 curated panels, the up- and down-regulated consensus components, "
     "and the 545 WGCNA intramodular hub genes."),
    ("T4", "Supplementary Table S4", "Table S4", ".xlsx, .csv",
     "Single-cell quality control and composition/state decomposition (GSE161277).",
     "Sub-table a: per-sample QC (raw and post-QC cell counts, epithelial fraction). "
     "Sub-table b: delta-log P, delta-log U and delta-log T for the Paneth signature per "
     "patient, with the decidable flag."),
    ("T5", "Supplementary Table S5", "Table S5", ".xlsx, .csv",
     "Four-cohort performance of the signature.",
     "Per-cohort AUC, same-direction concordance and platform annotation for each panel and "
     "scoring caliber; includes the GSE128435 exclusion record."),
    ("F1", "Fig. S1", "Figure S1", ".png, .pdf",
     "Global drift diagnostics.",
     "Genome-wide inter-group offsets per cohort and the residual offset remaining after "
     "sample-level correction."),
    ("F2", "Fig. S2", "Figure S2", ".png, .pdf",
     "Dual-caliber screening matrix (13 panels x 3 cohorts x 2 calibers).",
     "Pattern label (P / D / A / ?) of each panel in each cohort under the raw and the "
     "drift-corrected caliber."),
    ("F3", "Fig. S3", "Figure S3", ".png, .pdf",
     "GSE303070 (CosMx spatial atlas) cannot serve as an unbiased composition baseline.",
     "Field-of-view sampling composition and paired-donor epithelial fractions, showing "
     "selection of epithelium-rich fields."),
    ("T3", "Supplementary Table S3", "Table S3", ".xlsx, .csv",
     "Three-layer pattern matrix and sensitivity analysis.",
     "Sub-table a: pattern label of each panel across the three cohorts under both calibers. "
     "Sub-table b: DEFA5+DEFA6 and DEFA5/DEFA6/OLFM4/SPINK1 sensitivity results."),
]


def unify_additional_refs(body):
    """正文首见序 → Additional file N；返回 (改写后 body, 清单 HTML, 映射表)。

    BMC 硬规定：附加文件按正文引用顺序编号，正文以文件名显式引用。
    """
    # 三层模式表（原 Table S3）此前只在 Fig. 6 图注里被引用；BMC 要求"cited ... within
    # the main body of text"，故在三层模式句补一处正文引用（图号已在此前完成重排，故用 Fig. 6）。
    body = body.replace("(Fig. 6; Fig. S2):",
                        "(Fig. 6; Fig. S2; Supplementary Table S3):", 1)

    # 只在正文（Background → 附加文件清单之前）内统计首见序
    lo = body.find("Background</h2>")
    hi = body.find('Additional files</h2>')
    scope = body[lo:hi] if lo >= 0 and hi > lo else body

    order = []
    for key, marker, *_ in ADDITIONAL_ITEMS:
        pos = scope.find(marker)
        order.append((pos if pos >= 0 else 10 ** 9, key))
    order.sort()
    num = {k: i + 1 for i, (_, k) in enumerate(order)}

    # 先处理带子表后缀的形式，避免 "S4" 误吃 "S4b"
    body = body.replace("Supplementary Table S4a",
                        f"Additional file {num['T4']} (sub-table a)")
    body = body.replace("Supplementary Table S4b",
                        f"Additional file {num['T4']} (sub-table b)")
    # 通用替换
    for key, marker, *_ in ADDITIONAL_ITEMS:
        body = body.replace(marker, f"Additional file {num[key]}")
    body = body.replace("Supplementary Files", "additional files")
    body = body.replace("supplementary tables (S1-S5)", "the additional files")

    # 附图图注标题：Figure S1. → Additional file n (Figure S1).
    def _hdr(m):
        n = m.group(1)
        key = {"1": "F1", "2": "F2", "3": "F3"}[n]
        return f'<span class="fig">Additional file {num[key]} (Figure S{n}).</span>'
    body = re.sub(r'<span class="fig">Figure S([123])\.</span>', _hdr, body)

    ordered = sorted(
        (num[k], f"<li><b>Additional file {num[k]}</b>: {label} ({fmt}) &mdash; {title} "
                 f"{desc}</li>")
        for k, _, label, fmt, title, desc in ADDITIONAL_ITEMS)
    return body, "<ul>" + "".join(r for _, r in ordered) + "</ul>", num

# ---------------------------------------------------------------- 工具
def read(fn):
    return open(fn, encoding="utf-8").read()


def wc(s):
    """英文词数（与 BMC 口径接近：以字母开头的 token）"""
    return len(re.findall(r"[A-Za-z][A-Za-z\-']*", s))


# ---------------------------------------------------------------- 1. 引用转换
PART_RE = re.compile(
    r"^\s*([A-Z][A-Za-zÀ-ÿ\u00c0-\u017f\-']+?)"
    r"(?:\s+(?:et al\.|(?:&\s*amp;|&)\s+[A-Z][A-Za-zÀ-ÿ\-']+"
    r"|and\s+[A-Z][A-Za-zÀ-ÿ\-']+))?"
    r"\s*,\s*((?:19|20)\d{2}[a-z]?)\s*$")

def convert_parenthetical(h):
    """(McCorry et al., 2018; Kreis et al., 2024) → [1, 2]

    注意：HTML 实体 &amp; 自带分号，直接按 ';' 切分会把 "Tyler &amp; Tirosh, 2021"
    切成两段，故先用哨兵替换实体再切分。
    """
    stats = {"ok": 0, "skip": 0, "fail": [], "skipped_text": []}
    SENT = "\x00AMP\x00"
    # 括号里可能夹着 <span class="mono">GSE109308</span> 这类行内标签（数据集号 + 文献同括），
    # 纯文本分词会卡在标签上 → 先按"标签 + 半角分号"切，只丢标签、保留两侧分句。
    TAG = re.compile(r"</?span[^>]*>")

    def sub(m):
        inner = m.group(1)
        if not re.search(r"(?:19|20)\d{2}", inner):
            return m.group(0)
        norm = inner.replace("&amp;", SENT)
        norm = TAG.sub("", norm)
        parts = [p.strip().replace(SENT, "&amp;") for p in norm.split(";") if p.strip()]
        nums = []
        keeps = []   # 非引文分句（如数据集号）原样保留，避免被角标替换吞掉
        for p in parts:
            mm = PART_RE.match(p)
            if not mm:
                # 不含年份的分句（如裸数据集号）不是文献，保留原文
                if not re.search(r"(?:19|20)\d{2}", p):
                    keeps.append(p)
                    continue
                stats["skip"] += 1
                stats["skipped_text"].append(p)
                return m.group(0)
            key = (mm.group(1), mm.group(2))
            if key not in CITE_NUM:
                stats["fail"].append(p)
                return m.group(0)
            nums.append(CITE_NUM[key])
        if not nums:
            return m.group(0)
        stats["ok"] += 1
        # ★ 必须排序：编号是按**括号内原来的 author-year 书写顺序**收集的，
        #   而 author-year 顺序由作者姓氏字母序决定，与 Vancouver 的编号序无关。
        #   实测例子：正文写 "(Joo et al., 2009; Radeva et al., 2010; Pagnini et al., 2011;
        #   López-Arribillaga et al., 2021)"，对应编号 14/13/15/16 → 输出 [14, 13, 15, 16]，
        #   一个括号里的编号不升序，直接违反 Vancouver。
        #   去重同理：同一文献在同一个括号里出现两次一定是用例错误。
        cite = "[" + ", ".join(str(n) for n in sorted(set(nums))) + "]"
        # ★ 2026-10-20 踩坑：括号里常夹数据集号，写成
        #   (<span class="mono">GSE109308</span>; Bala et al., 2023)
        #   若整段替换成 [10]，GSE109308 会被**静默删掉**（构建不报错、自检也查不到，
        #   因为只查了角标而不查数据集号是否还在）。故保留所有非引文分句原文，只把
        #   引文分句换成角标。
        if not keeps:
            return cite
        return "(" + "; ".join(keeps + [cite]) + ")"

    out = re.sub(r"\(([^()]{3,200}?)\)", sub, h)
    return out, stats


def normalize_dataset_cites(h):
    """(<span class="mono">GSE109308</span>)<span class="small"> [10]</span>
       → <span class="mono">GSE109308</span> [10]

    ★ 2026-10-20 加。分节初稿的习惯写法是
        (<span class="mono">GSE109308</span>)(Bala et al., 2023)
      转换后变成 "(GSE109308) [10]" —— 数据集号被包在括号里、角标在括号外，
      这是 author-year 版式的残留，Vancouver 里读起来像"GSE109308 这个句子"。
      统一收成 "GSE109308 [10]"，与正文其他 GSE 号的书写方式一致。
    """
    pat = re.compile(
        r'\((<span class="mono">GSE\d+</span>)\)\s*<span class="small">\s*'
        r'(\[[\d,\s\-]+\])\s*</span>')
    n = len(pat.findall(h))
    h = pat.sub(r"\1 \2", h)
    return h, n


def convert_narrative(h):
    """X and colleagues → X et al. [n] ; Tyler and Tirosh → Tyler and Tirosh [3]"""
    stats = []

    def sub(m):
        name = m.group(1)
        if name in NARRATIVE:
            stats.append((name, NARRATIVE[name]))
            return f"{name} et al. [{NARRATIVE[name]}]"
        return m.group(0)

    h = re.sub(r"\b([A-Z][A-Za-zÀ-ÿ\-']+)\s+and\s+colleagues\b", sub, h)

    def sub2(m):
        stats.append(("Tyler and Tirosh", 3))
        return "Tyler and Tirosh [3]"

    h = re.sub(r"\bTyler and Tirosh\b(?!\s*\[)", sub2, h)

    # ---- 裸 (Author et al., YYYY) 兜底 ----
    # ★ 2026-10-20 发现的真 bug：分节初稿把一手文献引用写成
    #   <span class="small"> (Zheng et al., 2022)</span>
    #   convert_parenthetical 的正则 \(([^()]{3,200}?)\) 被 <span ...> 里的括号
    #   挡在标签之外 —— 括号仍能匹配，但 PART_RE 要求整段以姓氏开头，于是命中后
    #   落到"跳过"分支，结果是 5 处人名引用**静默留在投稿稿里**，而自检只查
    #   角标闭环（人名不算角标，查不出来）。
    #   兜底：整稿扫一遍裸形式 (X et al., YYYY) 与 (X & Y, YYYY)，能匹配的直接换角标。
    def sub3(m):
        key = (m.group(1), m.group(2))
        if key not in CITE_NUM:
            stats.append(("UNMAPPED", "%s %s" % (m.group(1), m.group(2))))
            return m.group(0)
        stats.append((m.group(1), CITE_NUM[key]))
        return "[%d]" % CITE_NUM[key]

    h = re.sub(
        r"\(([A-Z][A-Za-zÀ-ÿ\u00c0-\u017f\-']+)"
        r"(?:\s+(?:et\s+al\.|&\s+[A-Z][A-Za-zÀ-ÿ\-']+"
        r"|and\s+[A-Z][A-Za-zÀ-ÿ\-']+))?"
        r",\s*((?:19|20)\d{2})[a-z]?\)",
        sub3, h)
    return h, stats


def add_software_cites(h):
    """Methods 中的软件首见处补编号（WGCNA=17, Enrichr=18/19, MSigDB=20）"""
    hits = []
    pat = r"(weighted gene co-expression network analysis \(WGCNA\)) in R"
    if re.search(pat, h):
        h = re.sub(pat, r"\1 [17] in R", h, count=1)
        hits.append(17)
    s = "Functional enrichment was performed with Enrichr"
    if s in h:
        h = h.replace(s, s + " [18, 19]", 1)
        hits += [18, 19]
    if "MSigDB_Hallmark_2020" in h:
        h = h.replace("MSigDB_Hallmark_2020", "MSigDB_Hallmark_2020 [20]", 1)
        hits.append(20)
    return h, hits


# ---------------------------------------------------------------- 2. 图号重排
# 正文首见序：1 → 6 → 3 → 4 → 5 → 2 → 7 → 8   ⇒  仅 Fig. 2 ↔ Fig. 6 互换
def renumber_figures_in_text(h):
    h = h.replace("Fig. 6", "@@FIGX@@")
    h = h.replace("Fig. 2", "Fig. 6")
    h = h.replace("@@FIGX@@", "Fig. 2")
    h = h.replace("Figure 6.", "@@FIGY@@")
    h = h.replace("Figure 2.", "Figure 6.")
    h = h.replace("@@FIGY@@", "Figure 2.")
    return h


# ---------------------------------------------------------------- 内部代码清理
# 草稿期的 R1–R5（Results 分节码）与 M1–M12（Methods 分节码）不能出现在投稿稿里。
# 逐条给出等价改写，避免机械替换破坏语义。
INTERNAL_FIXES = [
    ("(see M3 and Results R1)", "(see Methods and Results)"),
    ("analysis (R2; Supplementary Table S1)", "analysis (Results; Supplementary Table S1)"),
    ("progenitor-like programme (R4), but", "progenitor-like programme (Results), but"),
    ("described in section R2.", "described above."),
    ("The principle underlying section R2 is established",
     "The principle underlying that analysis is established"),
    ("Our R2 finding", "Our own finding"),
    ("our single-cell design in R4 answers", "our single-cell design answers"),
    ("we adopt in R4 is likewise", "we adopt for the single-cell analysis is likewise"),
    ("The phenomenon quantified in R3 has", "The phenomenon quantified here has"),
    ("Our R3 finding that", "Our finding that"),
    ("what section R4 reports.", "what the single-cell analysis reports."),
    ("(see Results R1)", "(see Results)"),
    ("not epithelium-intrinsic (R2);", "not epithelium-intrinsic;"),
    ("survives all calibers (R3\u2013R4);", "survives all calibers;"),
    ("split into components (R5).", "split into components."),
    ("reporting components separately (R5).", "reporting components separately."),
    ("see Methods M10", "see Methods"),
    ("listed in M1", "listed in the Methods"),
    ("listed in the Methods (M1)", "listed in the Methods"),
    # Methods 里原有一份与 Declarations 重复的可用性声明；收口为指向 Declarations 的短句
    ("All datasets analysed here are publicly available from the NCBI Gene Expression "
     "Omnibus under the accessions listed in the Methods; no new data were generated. "
     "Analysis scripts, intermediate summary tables and the per-gene attribution table "
     "are available as Supplementary Files and, upon publication, from the corresponding "
     "repository. The signature gene list, the per-gene attribution table (313 genes "
     "&times; 15 columns), the panel-by-cohort discrimination matrix and the single-cell "
     "decomposition table are provided as Supplementary Tables S1&ndash;S5.",
     "All datasets analysed here are publicly available from the NCBI Gene Expression "
     "Omnibus (GEO) under the accessions listed in the Availability of data and materials "
     "section; no new data were generated. The derived tables and the analysis code are "
     "provided as additional files and in the code repository, as detailed in the "
     "Availability of data and materials section.")
]

# R5 两处前人工作论断原本无角标；经 PubMed 核实与已引文献逐字对应，补上编号
#   "spatial CODEX profiling of FAP mucosa" → Guha 2025 [8]
#     （其摘要原文："Our spatial CODEX experiments reveal that... FAP mucosa, pre-cancer
#       polyps and colorectal cancers exhibit substantial alterations in the cell type
#       composition and tissue microenvironment"）
#   "work in mouse and hereditary-polyposis models ... conserved events preceding
#    carcinoma development" → Bala 2023 [10]
#     （其摘要原文："Using genetically engineered and carcinogen-induced mouse models...
#       impaired differentiation is a conserved event preceding cancer development...
#       a patient with hereditary polyposis ... regenerative activity ... fetal genes"）
CITATION_COMPLETION = [
    ("conserved events preceding carcinoma development, and with spatial CODEX",
     "conserved events preceding carcinoma development [10], and with spatial CODEX"),
    ("alteration already in pre-cancerous polyps.",
     "alteration already in pre-cancerous polyps [8]."),
]


def clean_internal_codes(h):
    hits = []
    for a, b in INTERNAL_FIXES + CITATION_COMPLETION:
        if a in h:
            hits.append(a)
            h = h.replace(a, b)
    h, n = re.subn(r"(<h3>)M\d+\.\s*", r"\1", h)   # Methods 小节标题去 "M1. " 前缀
    return h, hits, n


# ---------------------------------------------------------------- 连续编号压缩
def compress_ranges(h):
    """Vancouver：连续编号用连字符区间（≥3 个连续才压缩），[1, 2, 3, 4] → [1-4]"""
    def sub(m):
        nums = [int(x) for x in re.findall(r"\d+", m.group(0))]
        out, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            if j - i >= 2:
                out.append(f"{nums[i]}-{nums[j]}")
            else:
                out.extend(str(x) for x in nums[i:j + 1])
            i = j + 1
        return "[" + ", ".join(out) + "]"
    return re.sub(r"\[\d+(?:\s*,\s*\d+)+\]", sub, h)


# ---------------------------------------------------------------- 主流程
def main():
    h = read(SRC)
    report = {}

    # ---- 切分章节 ----
    head_end = h.find('<h2 class="chap">')
    front = h[:head_end]
    CHAP_RE = re.compile(r'<h2 class="chap"><span class="chapno">(\d+)</span>(.*?)</h2>')
    ms = list(CHAP_RE.finditer(h[head_end:]))
    ms = [(m, m.start() + head_end, m.end() + head_end) for m in ms]
    chapters, tail = [], h[ms[-1][2]:]
    for i, (m, s, e) in enumerate(ms):
        nxt = ms[i + 1][1] if i + 1 < len(ms) else len(h)
        chapters.append(dict(num=int(m.group(1)),
                             title=re.sub(r"<[^>]+>", "", m.group(2)).strip(),
                             html=h[e:nxt]))
    report["chapters_in"] = [(c["num"], c["title"]) for c in chapters]

    # ---- 1. 清脚手架标签 ----
    n_lbl = 0
    def strip_scaffold(t):
        nonlocal n_lbl
        before = t
        t = re.sub(r'<span class="lbl">¶\d+</span>\s*—\s*', "", t)
        t = re.sub(r'<b>Paragraph \d+ — [^<]*</b>\s*<br>\s*', "", t)
        t = re.sub(r'<b>Paragraph \d+ — [^<]*</b>\s*', "", t)
        t = re.sub(r'<span class="tag"[^>]*>MS ¶\d+</span>\s*', "", t)
        if len(t) != len(before):
            n_lbl += 1
        return t
    for c in chapters:
        if c["title"] in ("Introduction", "Results", "Discussion"):
            c["html"] = strip_scaffold(c["html"])
    report["scaffold_labels_cleaned_chapters"] = n_lbl

    # ---- 2. 引用转换（全稿） ----
    for c in chapters:
        if c["title"] in ("Figure Legends", "References"):
            continue
        c["html"], st = convert_parenthetical(c["html"])
        report.setdefault("paren", {"ok": 0, "fail": [], "skip": 0, "skipped_text": []})
        report["paren"]["ok"] += st["ok"]
        report["paren"]["fail"] += st["fail"]
        report["paren"]["skip"] += st["skip"]
        report["paren"]["skipped_text"] += st["skipped_text"]
        c["html"], nar = convert_narrative(c["html"])
        report.setdefault("narrative", []).extend(nar)
        c["html"], n_ds = normalize_dataset_cites(c["html"])
        report["dataset_cites_unwrapped"] = report.get("dataset_cites_unwrapped", 0) + n_ds
    # 软件引用
    for c in chapters:
        if c["title"] == "Methods":
            c["html"], sw = add_software_cites(c["html"])
            report["software_cites"] = sw

    # ---- 2b. 连续编号压缩 + 内部代码清理 ----
    fixes_applied, h3_stripped = [], 0
    for c in chapters:
        c["html"] = compress_ranges(c["html"])
        c["html"], hits, n = clean_internal_codes(c["html"])
        fixes_applied += hits
        h3_stripped += n
    report["internal_fixes"] = fixes_applied
    report["methods_h3_prefix_stripped"] = h3_stripped

    # ---- 3. 图号重排 ----
    for c in chapters:
        c["html"] = renumber_figures_in_text(c["html"])

    # ---- 4. 组装 ----
    out = []

    # (a) Title page
    title = re.search(r"<h1>(.*?)</h1>", front, re.S).group(1).strip()
    out.append(f"<h1>{title}</h1>")
    out.append(title_page_html())
    out.append(TITLE_META_HTML)

    # (b) Abstract（BMC 三段式）
    total = sum(wc(v) for v in ABSTRACT.values())
    report["abstract_words"] = total
    out.append('<h2 class="chap"><span class="chapno">1</span>Abstract</h2>')
    out.append(f'<p class="meta">BMC structured abstract &mdash; '
               f'Background / Results / Conclusions &mdash; {total} words (limit 350) &mdash; '
               f'no literature citations.</p>')
    out.append('<div class="ms-para abstract">')
    for k in ("Background", "Results", "Conclusions"):
        out.append(f'<p><b class="abslbl">{k}</b> {ABSTRACT[k]}</p>')
    out.append("</div>")
    # Keywords
    out.append(f'<p class="kw"><b>Keywords:</b> {"<b>;</b> ".join(KEYWORDS)}</p>')

    # (c) 其余章节
    chap_map = {c["title"]: c for c in chapters}
    for c in chapters:
        if c["title"] == "Abstract":
            continue
        if c["title"] == "Introduction":
            out.append('<h2 class="chap"><span class="chapno">2</span>Background</h2>')
            out.append(c["html"].strip())
        elif c["title"] == "Methods":
            out.append('<h2 class="chap"><span class="chapno">3</span>Methods</h2>')
            out.append(c["html"].strip())
        elif c["title"] == "Results":
            out.append('<h2 class="chap"><span class="chapno">4</span>Results</h2>')
            out.append(add_results_headings(c["html"]))
        elif c["title"] == "Discussion":
            disc, concl = split_discussion(c["html"])
            out.append('<h2 class="chap"><span class="chapno">5</span>Discussion</h2>')
            out.append(disc)
            out.append('<h2 class="chap"><span class="chapno">6</span>Conclusions</h2>')
            out.append(concl)
        elif c["title"] == "Figure Legends":
            pass  # 稍后统一放（在 References 之后）
        elif c["title"] == "References":
            pass

    # (d) List of abbreviations
    out.append('<h2 class="chap"><span class="chapno">7</span>List of abbreviations</h2>')
    out.append('<div class="abbr"><ul>' + "".join(
        f"<li><b>{a}</b> &mdash; {b}</li>" for a, b in ABBREVIATIONS) + "</ul></div>")

    # (e) Declarations
    out.append('<h2 class="chap"><span class="chapno">8</span>Declarations</h2>')
    out.append('<div class="decl">')
    for k, v in DECLARATIONS:
        out.append(f"<h3>{k}</h3>")
        for para in v.split("\n"):
            out.append(f"<p>{para}</p>")
    out.append("</div>")

    # (f) References（编号 Vancouver）
    refs = json.load(open("_bmc_references.json", encoding="utf-8"))
    refs = sorted(refs, key=lambda r: r["num"])
    out.append('<h2 class="chap"><span class="chapno">9</span>References</h2>')
    out.append('<div class="sc-ref"><ol class="reflist">')
    for r in refs:
        out.append(f'<li value="{r["num"]}">{r["vancouver"]}</li>')
    out.append("</ol></div>")
    report["refs"] = len(refs)

    # (g) Figure legends（重编号 + 排序）
    legs, leg_order, leg_words = extract_legends(chap_map["Figure Legends"]["html"])
    out.append('<h2 class="chap"><span class="chapno">10</span>Figure legends</h2>')
    out.append('<div class="sc-leg">' + "\n".join(legs) + "</div>")
    report["legend_order"] = leg_order
    report["legend_words"] = leg_words

    # (h) Additional files —— 编号由正文首见序决定，故先占位、组装后再统一改写
    AF_SLOT = "@@ADDITIONAL_FILES_SLOT@@"
    out.append('<h2 class="chap"><span class="chapno">11</span>'
               'Additional files</h2>')
    out.append(AF_SLOT)

    body = "\n".join(out)
    body, af_list, af_num = unify_additional_refs(body)
    body = body.replace(AF_SLOT, f'<div class="addfiles">{af_list}</div>')
    report["additional_files"] = af_num
    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>{title[:90]}</title>
<style>{CSS}</style></head>
<body><div class="wrap">
{body}
<p class="small" style="margin-top:40px;border-top:1px solid var(--line);padding-top:14px">
BMC Genomics submission version &middot; generated 2026-10-16 by
<code>build_bmc_manuscript.py</code> &middot; English only.
</p>
</div></body></html>"""
    open(OUT, "w", encoding="utf-8").write(html)
    return html, report


# ---------------------------------------------------------------- Results 小节标题
R2_HEAD = ("Laser-capture microdissection assigns the tissue-level increase to stromal "
           "and myeloid composition")
R3_HEAD = "A three-layer pattern with a robust adenoma-stage Paneth-cell peak"
R5_HEAD = ("An independent spatial atlas cannot serve as a composition benchmark, and "
           "relation to previous work")


def add_results_headings(html):
    """在 R2 / R3 / R5 首段前插入小节标题。"""
    blocks = re.split(r'(?=<div class="ms-para">)', html)
    # 统计 ms-para 起点顺序
    idx = [i for i, b in enumerate(blocks) if b.startswith('<div class="ms-para">')]
    # 依据原 block 枚举：R1=15-19(5), R2=20-25(6), R3=26-30(5), R4=31-35(5), R5=36-41(6)
    r2_at, r3_at, r5_at = 5, 11, 21
    metas = [(r2_at, R2_HEAD), (r3_at, R3_HEAD), (r5_at, R5_HEAD)]
    for k, (pos, txt) in enumerate(sorted(metas, reverse=True)):
        blocks.insert(idx[pos], f"<h3>{txt}</h3>\n")
    return "".join(blocks)


# ---------------------------------------------------------------- Discussion 拆分
def split_discussion(html):
    blocks = re.split(r'(?=<div class="ms-para">)', html)
    idx = [i for i, b in enumerate(blocks) if b.startswith('<div class="ms-para">')]
    last = idx[-1]
    concl = blocks[last]
    disc = "".join(blocks[:last])
    return disc.strip(), concl.strip()


# ---------------------------------------------------------------- 图注提取
FIG_ORDER = {"Figure 1": 1, "Figure 2": 2, "Figure 3": 3, "Figure 4": 4,
             "Figure 5": 5, "Figure 6": 6, "Figure 7": 7, "Figure 8": 8,
             "Figure S1": 101, "Figure S2": 102, "Figure S3": 103}


def extract_legends(html):
    """图注块无嵌套 div，故首个 </div> 即为块边界。按新图号排序。"""
    idxs = [m.start() for m in re.finditer(r'<div class="leg">', html)] + [len(html)]
    blocks = []
    for i in range(len(idxs) - 1):
        b = html[idxs[i]:idxs[i + 1]]
        c = b.find("</div>")
        blocks.append(b[:c + 6].strip())

    def key(b):
        m = re.search(r"<h2>(.*?)</h2>", b, re.S)
        if not m:
            return 999
        name = re.sub(r"<[^>]+>", "", m.group(1)).strip().split(".")[0].strip()
        return FIG_ORDER.get(name, 998)

    blocks.sort(key=key)

    def headline(b):
        m = re.search(r"<h2>(.*?)</h2>", b, re.S)
        return re.sub(r"<[^>]+>", "", m.group(1)).strip() if m else "?"

    order = [headline(b)[:46] for b in blocks]
    words = {headline(b): wc(re.sub(r"<[^>]+>", " ", b)) for b in blocks}
    return blocks, order, words


# ---------------------------------------------------------------- CSS
CSS = """
:root{--blue:#2563eb;--teal:#0d9488;--amber:#d97706;--rose:#e11d48;
--violet:#7c3aed;--green:#16a34a;--ink:#1e293b;--mut:#64748b;--bg:#f8fafc;--line:#e2e8f0}
*{box-sizing:border-box}
body{font-family:"Segoe UI","Microsoft YaHei",sans-serif;background:var(--bg);color:var(--ink);
margin:0;padding:24px;line-height:1.75;font-size:15px}
.wrap{max-width:960px;margin:0 auto}
h1{font-size:23px;border-bottom:3px solid var(--blue);padding-bottom:10px;margin:0 0 10px;line-height:1.4}
h2{font-size:18px;color:var(--blue);margin:30px 0 10px}
h3{font-size:15.5px;color:var(--teal);margin:22px 0 8px}
p{margin:10px 0}
.meta{font-size:12.5px;color:var(--mut)}
.small{font-size:12.5px;color:var(--mut)}
.authline{font-size:14px;margin:4px 0}
h2.chap{font-size:20px;color:#fff;background:var(--ink);padding:11px 18px;border-radius:9px;
margin:44px 0 18px;border:0;line-height:1.45}
h2.chap .chapno{display:inline-block;background:rgba(255,255,255,.18);border-radius:6px;
padding:0 9px;margin-right:10px;font-size:16px;vertical-align:2px}
.ms-para{background:#fff;border:1px solid var(--line);border-radius:8px;padding:16px 20px;
margin:12px 0;font-family:Georgia,"Times New Roman",serif;font-size:15.5px;line-height:1.9}
.ms-para p{margin:12px 0}
.ms-para .num{color:var(--rose);font-weight:600}
.ms-para .mono,.mono{font-family:Consolas,monospace;font-size:13px;background:#f1f5f9;
padding:1px 6px;border-radius:4px}
.abstract{background:#f0fdf4;border-left:5px solid var(--green)}
.abslbl{color:var(--green);font-family:"Segoe UI",sans-serif;font-size:13.5px;
letter-spacing:.3px;margin-right:6px}
.kw{background:#eff6ff;border-left:4px solid var(--blue);border-radius:6px;padding:10px 14px;
font-size:13.5px}
.abbr,.addfiles{background:#fff;border:1px solid var(--line);border-radius:8px;padding:8px 18px}
.abbr li,.addfiles li{margin:3px 0;font-size:13.5px}
.decl h3{color:var(--ink);font-size:14.5px;margin:18px 0 4px;border-left:4px solid var(--teal);
padding-left:9px}
.decl p{font-size:14px;margin:4px 0 0;overflow-wrap:anywhere;word-break:break-word}
.sc-ref ol.reflist{padding-left:26px;font-size:13.5px}
.sc-ref ol.reflist li{margin:7px 0;line-height:1.6}
.sc-leg .leg{background:#fff;border:1px solid var(--line);border-radius:10px;padding:16px 20px;
margin:0 0 14px}
.sc-leg .leg h2{font-size:16px;margin:0 0 6px;color:var(--ink);border-left:4px solid var(--blue);
padding-left:9px}
.sc-leg .leg p{margin:6px 0}
.sc-leg .num,.sc-leg .k{color:var(--mut)}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.5px}
th{background:#1e293b;color:#fff;padding:8px 10px;text-align:left}
td{border:1px solid var(--line);padding:7px 10px;background:#fff;vertical-align:top}
tr:nth-child(even) td{background:#f8fafc}
@media print{
  body{background:#fff;padding:0;font-size:10.5pt}
  .wrap{max-width:none}
  .ms-para{font-size:10.5pt;line-height:1.45;padding:8px 12px;margin:8px 0;
    border:0;border-radius:0;background:#fff}
  .ms-para p{margin:7px 0}
  h1{font-size:15pt;margin-bottom:6px}
  h2.chap{font-size:12pt;padding:6px 12px;margin:22px 0 12px;border-radius:5px;
    break-after:avoid}
  h3{font-size:10.5pt;margin:14px 0 5px}
  .abbr,.decl,.addfiles,.kw{background:#fff;border:0;padding:0}
  .sc-leg .leg{padding:8px 12px;margin:0 0 8px}
  .sc-ref ol.reflist{font-size:10pt}
  .sc-ref ol.reflist li{margin:4px 0}
  .ms-para,.leg,.abbr,.decl,.addfiles{break-inside:auto}
  li,tr,blockquote{break-inside:avoid}
  h1,h2,h2.chap,h3{break-after:avoid}
  table{margin:8px 0}
  @page{margin:15mm 16mm}
}
"""

if __name__ == "__main__":
    html, rep = main()
    print("=== 章节（输入）", rep["chapters_in"])
    print(f"摘要词数          : {rep['abstract_words']}  (上限 350)")
    print(f"圆括号引用转换成功: {rep['paren']['ok']} 处；跳过 {rep['paren']['skip']} {rep['paren']['skipped_text']}；失败 {rep['paren']['fail']}")
    print(f"叙事式引用转换    : {len(rep['narrative'])} 处 -> {rep['narrative']}")
    print(f"软件引用插入      : {rep.get('software_cites')}")
    print(f"参考文献条数      : {rep['refs']}")
    print(f"清除脚手架标签章节: {rep['scaffold_labels_cleaned_chapters']}")
    print(f"内部代码改写      : {len(rep['internal_fixes'])} 条 -> {rep['internal_fixes']}")
    print(f"Methods 小节去码  : {rep['methods_h3_prefix_stripped']} 处")
    print(f"图注排序          : {rep['legend_order']}")
    print()
    print("图注词数：")
    for t, n in rep["legend_words"].items():
        print(f"   {n:4d}  {t}")

    # 残留内部代码自检
    body = re.sub(r"<style>.*?</style>", "", html, flags=re.S)
    body = body[:body.find("References</h2>")]
    resid_r = re.findall(r"\(?\bR[1-5]\b\)?", body)
    resid_m = re.findall(r"\bM(?:1[0-2]|[1-9])\b", body)
    print(f"\n残留 R 码: {resid_r if resid_r else '无'}")
    print(f"残留 M 码: {resid_m if resid_m else '无'}")
    print(f"\n写出 {OUT}  ({os.path.getsize(OUT):,} bytes)")
