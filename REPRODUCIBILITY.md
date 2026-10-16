# Reproducibility manifest — adenoma–carcinoma compartment attribution study

This file is the supporting material for the **Availability of data and materials** and
**Code availability** statements of the manuscript submitted to *BMC Genomics*.
It documents the full chain: public input data → analysis scripts → derived tables →
figures → manuscript.

> **Status of the mapping.** Every `script → output` pair below was verified by reading the
> script's own input/output filenames at the time of release. The *order* of the stages is a
> logical reconstruction of the pipeline; **[Author: confirm the execution order before
> archiving]**.

---

## 1. Environment

| Item | Value |
|---|---|
| Operating system | Windows 11 (code is platform independent) |
| Python | 3.13 |
| Python packages | `numpy`, `pandas`, `scipy`, `openpyxl`, `requests`, `matplotlib` |
| R | `limma` (differential expression), `WGCNA` (co-expression network) — [Author: confirm R version] |
| External services | NCBI E-utilities (literature), Enrichr, MSigDB (`MSigDB_Hallmark_2020`) |
| Hardware notes | The 10x single-cell objects and the CosMx object were processed on a 16 GB workstation |

A frozen environment file (`REQUIREMENTS.txt`) is provided at the repository root, so the
environment can be rebuilt from scratch.

---

## 2. Public input data (no new data were generated)

All raw data are pre-existing and public. Accessions, with permanent identifier:

| GEO accession | Role in the study | URL |
|---|---|---|
| GSE41258 | Discovery cohort (whole tissue, Affymetrix) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258 |
| GSE37364 | Discovery cohort (whole tissue) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37364 |
| GSE20916 | Discovery cohort (whole tissue) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20916 |
| GSE41657 | Discovery cohort (whole tissue) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41657 |
| GSE4183 | Three-layer validation (whole tissue, Affymetrix) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4183 |
| GSE117606 | Three-layer validation (whole tissue, Agilent) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117606 |
| GSE15960 | LCM pure-epithelium cohort (Affymetrix) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15960 |
| GSE128435 | Screened then **excluded** (composition-confounded); recorded for completeness | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128435 |
| GSE161277 | Single-cell cohort (normal → adenoma → carcinoma) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161277 |
| GSE303070 | CosMx spatial atlas — used only as a descriptive counter-example | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303070 |

All URLs accessed 16 October 2026. Acquisition helpers: `geo_dl.py`, `geo_dl2.py`,
`dl_ebi.py` / `dl_loop.py` / `dl_persist.py` (ArrayExpress mirrors), `sc161277_dl_final.py`
(single-cell), `gse303070_probe.py` / `gse303070_tissue.py` (spatial).

---

## 3. Pipeline stages

### Stage A — consensus signature (313 genes)
```
consensus_signature.py   ->  consensus_all.csv, consensus_core.csv, consensus_summary.json
sig_validate.py          ->  cross-platform validation of the signature
```
Stouffer integration of per-cohort differential calls.
Reads the three per-cohort DEG tables in `02.result/` plus three files from `upstream_r/`
(§5): `01.limma/DEG_all_Tumor_vs_Polyp.csv`, `02.WGCNA/WGCNA_module_genes.csv`,
`03.enrich/gene_Tumor_vs_Polyp.csv`.

### Stage B — drift diagnostics
```
drift_diag.py            ->  drift_diag.json      # genome-wide inter-group offsets per cohort
lcm_diag.py              ->  LCM cohort offset diagnostics
```
Output feeds main text (dual-caliber reporting) and **Additional file 5 (Figure S1)**.

### Stage C — three-layer panel screening
```
three_layer_check.py     ->  three_layer_pattern.json, three_layer_results.json
r3_sens_paneth.py        ->  r3_sens_paneth.json   # DEFA5+DEFA6 / 4-gene sensitivity
r3_facts.py              ->  r3_facts.json
```
Feeds **main Figures 5–6** and **Additional file 7 (Table S3)** + **Additional file 1 (Table S1)**.

### Stage D — single-cell decomposition (T = P × U)
```
sc161277_scan.py         ->  sc161277_scan.json
sc161277_pairs.py        ->  sc161277_pairs.json   # per-patient P / U / T terms
```
Feeds **main Figures 7–8** and **Additional file 3 (Table S4)**.

### Stage E — spatial reference (descriptive only)
```
sc303070_compo.py / _compo2.py / _compo3.py   ->  sc303070_composition_final.json
sc303070_niche.py                             ->  sc303070_niche.json
```
Feeds **Additional file 8 (Figure S3)** — documents the field-of-view selection bias.

### Stage F — gene-set membership and enrichment
```
enrichr_dual.py          ->  enrichr_results.json     # Enrichr (GO / MSigDB Hallmark)
gse117606_paired.py      ->  gse117606_paired_results.json
panels_def.json          ->  panel definitions (17 curated panels)
```
Feeds **Additional file 2 (Table S2)**.

### Stage G — derived tables
```
build_suppl_tables.py
    reads : suppl_table1_313genes.csv, panels_def.json, r3_facts.json,
            r3_sens_paneth.json, sc161277_scan.json, sc161277_pairs.json,
            manuscript_facts.json,
            upstream_r/02.WGCNA/WGCNA_top_module_genes.csv    (see §5)
            upstream_r/03.enrich/Hub_Tumor_vs_Polyp.csv       (see §5)
    writes: suppl/S1_313gene_attribution.csv
            suppl/S2_panel_gene_membership.csv
            suppl/S3_panel_pattern_matrix.csv
            suppl/S3b_paneth_sensitivity.csv
            suppl/S4a_sc161277_sample_QC.csv
            suppl/S4b_sc161277_decomposition.csv
            suppl/S5_four_cohort_performance.csv
            suppl/Supplementary_Tables_S1-S5.xlsx
```

### Stage H — figures
```
make_figs.py        ->  fig1..fig6 (study design, three-layer pattern, 2x2 decomposition,
                        per-gene scatter, panel x cohort AUC heatmap, matched-pair bars)
make_figs_sc.py     ->  fig7_单细胞P×U分解.png, fig8_DEFA腺瘤峰dot.png
make_figS1.py       ->  figS1_漂移诊断.png
make_figS3.py       ->  figS3_GSE303070取样偏差.png
```

### Stage I — manuscript build and QA
```
build_bmc_manuscript.py   ->  Manuscript_BMC_Genomics_投稿版.html   (BMC Genomics format)
_fetch_bmc_refs.py        ->  _bmc_references.json                   (Vancouver, re-fetched authors)
verify_bmc.py             ->  19-point compliance self-check (all pass)
gen_bmc_readiness.py      ->  BMC投稿准备报告_BMC_Genomics.html
```

---

## 4. Reproduce one artefact end-to-end (example)

```bash
# Additional file 5 (Figure S1) — global drift diagnostics
python drift_diag.py           # -> drift_diag.json
python make_figS1.py           # -> figs/figS1_漂移诊断.png

# Additional file 7 (Table S3) — three-layer pattern matrix
python three_layer_check.py    # -> three_layer_pattern.json, three_layer_results.json
python r3_sens_paneth.py       # -> r3_sens_paneth.json
python build_suppl_tables.py   # -> suppl/S3_panel_pattern_matrix.csv, S3b_paneth_sensitivity.csv

# Main Figures 7-8 — single-cell decomposition
python sc161277_scan.py        # -> sc161277_scan.json
python sc161277_pairs.py       # -> sc161277_pairs.json
python make_figs_sc.py         # -> figs/fig7_单细胞P×U分解.png, figs/fig8_DEFA腺瘤峰dot.png
```

---

## 5. Upstream R prerequisites (shipped in `upstream_r/`)

Two scripts consume CSV files produced by the upstream **R** pipeline (`limma` differential
expression, `WGCNA` co-expression network, enrichment intersection). These five files are
**shipped with the repository** under `upstream_r/`, so the Python package is self-contained:

| Consumed by | Original path | Shipped as |
|---|---|---|
| `consensus_signature.py` | `D:/2/01.limma/DEG_all_Tumor_vs_Polyp.csv` | `upstream_r/01.limma/DEG_all_Tumor_vs_Polyp.csv` |
| `consensus_signature.py` | `D:/2/02.WGCNA/WGCNA_module_genes.csv` | `upstream_r/02.WGCNA/WGCNA_module_genes.csv` |
| `consensus_signature.py`, `build_suppl_tables.py` | `D:/2/02.WGCNA/WGCNA_top_module_genes.csv` | `upstream_r/02.WGCNA/WGCNA_top_module_genes.csv` |
| `consensus_signature.py` | `D:/2/03.enrich/gene_Tumor_vs_Polyp.csv` | `upstream_r/03.enrich/gene_Tumor_vs_Polyp.csv` |
| `build_suppl_tables.py` | `D:/2/03.enrich/Hub_Tumor_vs_Polyp.csv` | `upstream_r/03.enrich/Hub_Tumor_vs_Polyp.csv` |

Both scripts resolve the location through the `AACA_UPSTREAM_R` environment variable,
defaulting to `./upstream_r` and falling back to the original `D:/2` path on the developer's
machine. Without these five files, **Additional file 2 cannot be rebuilt** — the chain breaks
at Stage A / Stage G.

> **Verified 16 October 2026.** Rebuilding `suppl/` from `upstream_r/` inside the repository
> reproduces the released tables exactly (all seven CSV md5 hashes unchanged).

---

## 6. Deposit checklist (BMC requirement)

BMC requires the code to be archived in a repository with a **DOI** (Zenodo is recommended
for GitHub projects).

- [x] Public repository assembled: this file, the Stage A–I scripts, `panels_def.json`,
      `suppl_table1_313genes.csv`, `manuscript_facts.json`, `r3_facts.json`,
      `r3_sens_paneth.json`, `three_layer_pattern.json`, `sc161277_scan.json`,
      `sc161277_pairs.json`, `drift_diag.json`, `sc303070_niche.json`,
      `gse117606_paired_results.json`, `enrichr_results.json`, the three per-cohort DEG
      tables under `02.result/`, and the five upstream R CSVs under `upstream_r/`.
- [x] Frozen environment recorded and committed (`REQUIREMENTS.txt`).
- [x] License chosen and stated (`LICENSE`, MIT).
- [x] Package verified self-contained: rebuilt from `upstream_r/` with identical CSV hashes.
- [ ] Push the repository and archive a tagged release on **Zenodo** to obtain a DOI.
      Automated by `push_github.py` (push → tag `v1.0.0` → GitHub release → retrieve the
      Zenodo DOI → write `code_availability.json` → rebuild the manuscript). The
      Zenodo–GitHub switch for this repository must be turned **ON before** the release is
      published, otherwise that release is not archived.
- [ ] Repository URL and Zenodo DOI into the *Code availability* block of the manuscript.
      Not done by hand: `build_bmc_manuscript.py` reads them from `code_availability.json`,
      which `push_github.py` writes.
- [ ] The author metadata in `.zenodo.json` (creators / affiliation) must match the
      manuscript title page — `push_github.py` refuses to publish the release while the
      `__FILL_ME__` placeholder is still present.
- [ ] Confirm the R version used for `limma` / `WGCNA` (`push_github.py --r-version "R x.y.z"`).

---

## 7. What is deliberately *not* reproducible from code alone

- The **GSE128435 exclusion** decision is a documented judgement call (composition
  confounding), recorded in Additional file 4; it is not an automatic filter.
- The **"non-arbitrable reference"** status of GSE303070 is a reasoned decision based on
  field-of-view selection, not a numeric threshold.
- All p-values are reported as **agreement between two exchangeable calibers** rather than as
  single significance tests; the criterion is stated in the Methods.

---

## 8. Byte-level determinism — what to hash and what not to

Checked empirically on the developer's machine:

| Output | Byte-identical across runs? | Reason |
|---|---|---|
| `suppl/*.csv` | **yes** | plain text, deterministic ordering |
| `figs/*.png` | **yes, if the Matplotlib version matches** | the only text chunk is `Software: Matplotlib version X.Y.Z`; no timestamp is embedded |
| `suppl/Supplementary_Tables_S1-S5.xlsx` | **no** | `openpyxl` writes `docProps/core.xml` with `dcterms:created` / `dcterms:modified`, and the ZIP entry headers carry the run time. Every cell is identical; only the container bytes differ. Two consecutive runs gave different MD5s. |

**Consequence for anyone verifying this deposit:** compare the **`.csv` tables**, not the
`.xlsx` workbook or a single archive checksum. A differing `.xlsx` hash does **not** indicate
a different result — the `.csv` files are authoritative and the `.xlsx` is a
reader-friendly repackaging of exactly those tables.

