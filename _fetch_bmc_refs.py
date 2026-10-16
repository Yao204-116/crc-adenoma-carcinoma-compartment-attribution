# -*- coding: utf-8 -*-
"""
按 BMC Genomics（BioMed Central Vancouver）参考文献格式重排 21 条引文。

BMC 官方示例规则（逐字）：
  * "List authors in the order they appear in the original publication."
  * "If there are more than six authors, list the first six authors followed by et al."
  * "Use the official PubMed journal title abbreviation."
  * "Journal Abbreviation. Year;Volume(Issue):Page range."
  * 期刊文章不带 DOI（官方示例均无 DOI）

引用编号 = 正文首次出现顺序：
  1 McCorry 2018 / 2 Kreis 2024 / 3 Tyler 2021 / 4 Block 2019 / 5 Yan 2022 /
  6 Becker 2022 / 7 Zheng 2022 / 8 Guha 2025 / 9 Yang 2025 / 10 Bala 2023 /
  11 Porter 1997 / 12 Mallow 1996 / 13 Radeva 2010 / 14 Joo 2009 / 15 Pagnini 2011 /
  16 López-Arribillaga 2021 / 17 Langfelder 2008 (WGCNA) / 18 Kuleshov 2016 (Enrichr) /
  19 Chen 2013 (Enrichr) / 20 Liberzon 2015 (MSigDB) / 21 Thiruvengadam 2019
"""
import json, re, subprocess, sys, time
import xml.etree.ElementTree as ET

ORDER = [
    (1,  "30105762", "McCorry",           "McCorry et al., 2018"),
    (2,  "38349551", "Kreis",             "Kreis et al., 2024"),
    (3,  "33972543", "Tyler",             "Tyler & Tirosh, 2021"),
    (4,  "31780722", "Block",             "Block et al., 2019"),
    (5,  "35433435", "Yan",               "Yan et al., 2022"),
    (6,  "35726067", "Becker",            "Becker et al., 2022"),
    (7,  "35221332", "Zheng",             "Zheng et al., 2022"),
    (8,  "39605357", "Guha",              "Guha et al., 2025"),
    (9,  "40075322", "Yang",              "Yang et al., 2025"),
    (10, "36989360", "Bala",              "Bala et al., 2023"),
    (11, "9169779",  "Porter",            "Porter et al., 1997"),
    (12, "8626737",  "Mallow",            "Mallow et al., 1996"),
    (13, "20979654", "Radeva",            "Radeva et al., 2010"),
    (14, "19269007", "Joo",               "Joo et al., 2009"),
    (15, "21346603", "Pagnini",           "Pagnini et al., 2011"),
    (16, "34831152", "López-Arribillaga", "López-Arribillaga et al., 2021"),
    (17, "19114008", "Langfelder",        "WGCNA (Langfelder & Horvath, 2008)"),
    (18, "27141961", "Kuleshov",          "Enrichr 2016 (Kuleshov et al., 2016)"),
    (19, "23586463", "Chen",              "Enrichr 2013 (Chen et al., 2013)"),
    (20, "26771021", "Liberzon",          "MSigDB Hallmark (Liberzon et al., 2015)"),
    (21, "31211760", "Thiruvengadam",     "Thiruvengadam et al., 2019"),
]

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"


def fetch(pmids):
    url = f"{BASE}/efetch.fcgi?db=pubmed&id={','.join(pmids)}&retmode=xml"
    out = subprocess.run(
        ["curl", "-sS", "--ssl-no-revoke", "--max-time", "90", "-A", UA, url],
        capture_output=True)
    if out.returncode != 0:
        sys.exit("curl failed: " + out.stderr.decode("utf-8", "replace"))
    return out.stdout.decode("utf-8", "replace")


def normalize(s):
    """统一异体撇号/连字符，避免投稿系统乱码。"""
    if not s:
        return ""
    for a, b in [("\u02bc", "'"), ("\u2019", "'"), ("\u2018", "'"),
                 ("\u2013", "-"), ("\u2014", "-"), ("&amp;", "&")]:
        s = s.replace(a, b)
    return re.sub(r"\s+", " ", s).strip()


def abbreviate_pages(p):
    """NLM/Vancouver 页码缩写：仅当首尾页码位数相同时才缩写，否则原样保留。

    422-426 -> 422-6 ; 2389-2395 -> 2389-95 ; 872-880 -> 872-80 ;
    4038-45 (已是缩写) -> 保持 4038-45
    """
    if not p:
        return ""
    p = p.replace("–", "-").strip()
    m = re.match(r"^([A-Za-z]*)(\d+)\s*-\s*([A-Za-z]*)(\d+)$", p)
    if not m:
        return p
    pre_a, a, pre_b, b = m.groups()
    if pre_a != pre_b or len(a) != len(b):
        return p
    i = 0
    while i < min(len(a), len(b)) and a[i] == b[i]:
        i += 1
    tail = b[i:]
    if not tail:
        return f"{a}-{b}"
    return f"{a}-{tail}"


def text(el, path, default=""):
    e = el.find(path)
    if e is None:
        return default
    return "".join(e.itertext()).strip()


def parse(xml):
    root = ET.fromstring(xml)
    out = {}
    for art in root.iter("PubmedArticle"):
        pmid = text(art, ".//MedlineCitation/PMID")
        a = art.find(".//MedlineCitation/Article")
        if a is None:
            continue
        journal = a.find("Journal")
        ta = ""
        if journal is not None:
            ta = (text(journal, "MedlineTA") or text(journal, "ISOAbbreviation")
                  or text(journal, "Title"))
        ta = normalize(ta)
        ji = a.find("Journal/JournalIssue")
        vol = text(ji, "Volume") if ji is not None else ""
        iss = text(ji, "Issue") if ji is not None else ""
        year = ""
        if ji is not None:
            year = text(ji, "PubDate/Year") or text(ji, "PubDate/MedlineDate")
        if not year:
            year = text(art, ".//ArticleDate/Year") or text(art, ".//PubMedPubDate/Year")
        year = re.sub(r"[^0-9].*$", "", year)
        title = normalize(re.sub(r"\s+", " ", text(a, "ArticleTitle")))
        title = re.sub(r"\.$", "", title)
        pages = text(a, "Pagination/MedlinePgn") or text(a, "Pagination/StartPage")
        pages = normalize(pages)
        artnum = ""
        for aid in a.iter("ArticleId"):
            if aid.get("IdType") in ("pii",) and aid.text:
                artnum = normalize(aid.text)
        for eid in a.iter("ELocationID"):
            if eid.get("EIdType") == "pii" and eid.text and not artnum:
                artnum = normalize(eid.text)
        doi = ""
        for eid in a.iter("ELocationID"):
            if eid.get("EIdType") == "doi" and eid.text:
                doi = normalize(eid.text)
        authors = []
        for au in a.iter("Author"):
            ln = text(au, "LastName")
            ini = text(au, "Initials")
            coll = text(au, "CollectiveName")
            if coll:
                authors.append(normalize(coll))
            elif ln:
                authors.append(normalize(f"{ln} {ini}"))
        # e-pub date fallback for year
        out[pmid] = dict(pmid=pmid, journal_ta=ta, volume=vol, issue=iss, year=year,
                         title=title, pages=pages, artnum=artnum, doi=doi,
                         authors=authors, n_authors=len(authors))
    return out


def fmt(rec):
    """BMC / Vancouver 期刊文章格式。

    有卷期页 -> "Journal. Year;Volume(Issue):Pages."
    无卷期页（如预印本） -> "Journal. Year. doi:..."; BMC 官方示例 "Article within a journal by DOI"
    """
    au = rec["authors"]
    if len(au) > 6:
        astr = ", ".join(au[:6]) + ", et al."
    else:
        astr = ", ".join(au)
    if not astr.endswith("."):
        astr += "."
    j = rec["journal_ta"]
    y, v, i, p = rec["year"], rec["volume"], rec["issue"], rec["pages"]
    if v:
        loc = f"{y};{v}"
        if i:
            loc += f"({i})"
        pg = abbreviate_pages(p) or rec["artnum"]
        if pg:
            loc += f":{pg}"
        loc += "."
    else:
        loc = f"{y}."
        if rec["doi"]:
            loc += f" doi:{rec['doi']}."
    s = f"{astr} {rec['title']}. {j}. {loc}"
    return re.sub(r"\s+", " ", s).strip()


if __name__ == "__main__":
    pmids = [p for _, p, _, _ in ORDER]
    xml = fetch(pmids)
    open("_bmc_refs_raw.xml", "w", encoding="utf-8").write(xml)
    recs = parse(xml)
    print(f"fetched {len(recs)}/{len(pmids)} records\n")
    result = []
    for num, pmid, _surname, intext in ORDER:
        r = recs.get(pmid)
        if not r:
            print(f"  !! MISSING {num} {pmid} {intext}")
            continue
        r["num"] = num
        r["intext"] = intext
        r["vancouver"] = fmt(r)
        result.append(r)
        print(f"{num:2d}. [PMID {pmid}] n_authors={r['n_authors']:2d}")
        print(f"    {r['vancouver']}")
    json.dump(result, open("_bmc_references.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"\nwrote _bmc_references.json ({len(result)} entries)")
