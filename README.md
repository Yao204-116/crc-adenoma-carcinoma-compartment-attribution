# Composition or cell state? — reproducibility package

Code and derived data supporting the manuscript submitted to *BMC Genomics*:

> **Composition or cell state? A three-resolution decomposition of the adenoma–carcinoma
> transcriptome identifies an epithelial-intrinsic adenoma-stage differentiation peak**

This repository contains **no new primary data**. Every script operates on public datasets
from NCBI GEO and on the derived tables that the earlier stages write to disk. The goal is
that each number, table and figure in the manuscript can be traced back to a script and an
input file.

---

## 1. Data sources

All expression data are pre-existing and public. No new specimens were collected.

| GEO accession | Role |
|---|---|
| GSE41258 | Discovery cohort (whole tissue, Affymetrix) |
| GSE37364 | Discovery cohort (whole tissue) |
| GSE20916 | Discovery cohort (whole tissue, macrodissected) |
| GSE41657 | Discovery cohort (whole tissue, Agilent) |
| GSE4183 | Three-layer validation (whole tissue, Affymetrix) |
| GSE117606 | Three-layer validation (whole tissue, Agilent) |
| GSE15960 | LCM pure-epithelium cohort (Affymetrix) |
| GSE161277 | Single-cell cohort (normal → adenoma → carcinoma) |
| GSE128435 | Screened, then **excluded** (composition-confounded); kept for transparency |
| GSE303070 | CosMx spatial atlas — used only as a descriptive counter-example |

Download helpers: `geo_dl.py`, `geo_dl2.py`, `dl_ebi.py` / `dl_loop.py` / `dl_persist.py`
(ArrayExpress mirrors), `sc161277_dl_final.py` (single-cell), `gse303070_probe.py` /
`gse303070_tissue.py` (spatial).

---

## 2. Environment

| Item | Value |
|---|---|
| Python | 3.9+ (developed on 3.13) |
| Python packages | see `REQUIREMENTS.txt` (`numpy`, `pandas`, `scipy`, `matplotlib`, `openpyxl`, `requests`) |
| R | `limma` and `WGCNA` — used **upstream** for the modules shipped in `upstream_r/` (see §5) |
| External services | NCBI E-utilities, Enrichr, MSigDB (`MSigDB_Hallmark_2020`) |

```bash
pip install -r REQUIREMENTS.txt
```

`geo_cache/` and platform annotation files (`GPL*.annot.gz`, `Homo_sapiens.gene_info.gz`)
are **not** shipped: the download scripts fetch them on first run and cache them locally.

---

## 3. Layout

```
.
├── *.py                    analysis scripts (Stages A–I, flat)
├── *.json / *.csv          derived intermediate results (the numbers behind the figures)
├── 02.result/              per-cohort DEG tables for GSE37364 / GSE20916 / GSE41657
├── upstream_r/             outputs of the upstream R pipeline required to rebuild tables
│   ├── 01.limma/           DEG_all_Tumor_vs_Polyp.csv
│   ├── 02.WGCNA/           WGCNA_module_genes.csv, WGCNA_top_module_genes.csv
│   └── 03.enrich/          gene_Tumor_vs_Polyp.csv, Hub_Tumor_vs_Polyp.csv
├── figs/                   rendered main and supplementary figures
├── suppl/                  Supplementary Tables S1–S5 (.csv and .xlsx)
└── REPRODUCIBILITY.md      the full script → output → manuscript-artefact mapping
```

Scripts locate their inputs relative to their own directory, so the package runs from any
working directory. Two optional environment variables override this:

| Variable | Effect |
|---|---|
| `AACA_BASE` | Root directory of the scripts (default: the script's own directory) |
| `AACA_UPSTREAM_R` | Location of the upstream R products (default: `./upstream_r`, falling back to `D:/2` on the developer's machine) |

---

## 4. Quick start — rebuild one artefact end to end

```bash
# Additional file 5 (Figure S1) — global drift diagnostics
python drift_diag.py            # -> drift_diag.json
python make_figS1.py            # -> figs/figS1_漂移诊断.png

# Main Figures 7–8 — single-cell decomposition (T = P x U)
python sc161277_scan.py         # -> sc161277_scan.json
python sc161277_pairs.py        # -> sc161277_pairs.json
python make_figs_sc.py          # -> figs/fig7_*.png, figs/fig8_*.png

# Supplementary Tables S1–S5
python build_suppl_tables.py    # -> suppl/*.csv, suppl/Supplementary_Tables_S1-S5.xlsx
```

`build_suppl_tables.py` needs `openpyxl` for the `.xlsx` output.

---

## 5. Upstream dependency (important)

`consensus_signature.py` and `build_suppl_tables.py` consume five CSV files produced by an
upstream **R** pipeline (`limma` differential expression, `WGCNA` co-expression network,
enrichment intersection). Those five files are shipped in `upstream_r/` so that the Python
package is self-contained; the R scripts themselves are not part of this repository.

Without them, **Additional file 2 cannot be rebuilt** — the chain breaks at Stage A / G.

---

## 6. What is deliberately *not* reproducible from code alone

- The **GSE128435 exclusion** is a documented judgement call (composition confounding), not
  an automatic filter.
- The **non-arbitrable reference** status of GSE303070 is a reasoned decision based on
  field-of-view selection bias, not a numeric threshold.
- All p-values are reported as **agreement between two exchangeable calibers** rather than as
  single significance tests; the criterion is stated in the Methods of the manuscript.

See `REPRODUCIBILITY.md` for the complete stage-by-stage mapping.

---

## 7. Citation, archiving and license

If you use this code, please cite the manuscript above (citation details will be updated
once the article is published).

| | |
|---|---|
| Source repository | \<REPO_URL\> |
| Archived version (DOI) | minted by the Zenodo–GitHub integration for every tagged release |
| License | **MIT** — see `LICENSE` |

Each tagged release of this repository is automatically archived by Zenodo, which mints a
citable DOI. The DOI of release `v1.0.0` is the one quoted in the *Code availability*
statement of the manuscript. Metadata for the DOI record is declared in `.zenodo.json`.

### Reproducing the archived release

```bash
git clone <REPO_URL>
cd <repo>
pip install -r REQUIREMENTS.txt
python build_suppl_tables.py      # smoke test: rebuilds Supplementary Tables S1-S5
```

