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
fig_en.py           shared English label map + render-time translation, and the
                    save-time guards (CJK residue, horizontal overflow, label
                    overlap).  Imported by all four scripts below.
check_fig_english.py static check: no on-screen Chinese left in the plotting scripts
make_figs.py        ->  fig2_56_pair_concordance.png,
                        fig3_2x2_decomposition_whole_tissue_vs_LCM.png,
                        fig4_per_gene_attribution_scatter.png,
                        fig5_panel_cohort_AUC_heatmap.png,
                        fig6_three_layer_pattern_adenoma_peak.png,
                        figS2_double_caliber_screening_matrix.png
make_figs_sc.py     ->  fig7_single_cell_PxU_decomposition.png,
                        fig8_DEFA6_adenoma_peak_dot.png
make_figS1.py       ->  figS1_drift_diagnostics.png
make_figS3.py       ->  figS3_GSE303070_sampling_bias.png
```
`figs/fig1_study_design.png` is rendered from the hand-authored
`figs/fig1_study_design.html` (a study-design schematic, not data-driven), so it has
no generator script.

Figure numbers follow first-mention order in the manuscript.  The file names under
`figs/` already carry that final number, so `figs/figN_*.png` maps one-to-one onto
`Figure_N.png` in the upload package (`make_upload_package.py`).

### Stage I — manuscript build and QA
```
build_bmc_manuscript.py   ->  Manuscript_BMC_Genomics_投稿版.html   (BMC Genomics format)
_fetch_bmc_refs.py        ->  _bmc_references.json                   (Vancouver, re-fetched authors)
verify_bmc.py             ->  21-point compliance self-check (all pass)
gen_bmc_readiness.py      ->  BMC投稿准备报告_BMC_Genomics.html
```

---

## 4. Reproduce one artefact end-to-end (example)

```bash
# Additional file 5 (Figure S1) — global drift diagnostics
python drift_diag.py           # -> drift_diag.json
python make_figS1.py           # -> figs/figS1_drift_diagnostics.png

# Additional file 7 (Table S3) — three-layer pattern matrix
python three_layer_check.py    # -> three_layer_pattern.json, three_layer_results.json
python r3_sens_paneth.py       # -> r3_sens_paneth.json
python build_suppl_tables.py   # -> suppl/S3_panel_pattern_matrix.csv, S3b_paneth_sensitivity.csv

# Main Figures 7-8 — single-cell decomposition
python sc161277_scan.py        # -> sc161277_scan.json
python sc161277_pairs.py       # -> sc161277_pairs.json
python make_figs_sc.py         # -> figs/fig7_single_cell_PxU_decomposition.png,
                               #    figs/fig8_DEFA6_adenoma_peak_dot.png
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
      Automated by `push_github.py` (push → tag `v1.0.2` → GitHub release → retrieve the
      Zenodo DOI → write `code_availability.json` → rebuild the manuscript). The
      Zenodo–GitHub switch for this repository must be turned **ON before** the release is
      published, otherwise that release is not archived.
- [ ] Repository URL and Zenodo DOI into the *Code availability* block of the manuscript.
      Not done by hand: `build_bmc_manuscript.py` reads them from `code_availability.json`,
      which `push_github.py` writes.
- [ ] The author metadata in `.zenodo.json` (creators / affiliation) must match the
      manuscript title page — `preflight_release.py` and `push_github.py` both refuse to
      publish the release while the creator placeholder is still present.
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
| `figs/*.png` | **only on the same binary stack** | layout depends on the FreeType build Matplotlib links against: the same Matplotlib 3.10.6 gave FreeType 2.13.3 under conda (matching) and 2.6.1 from the PyPI wheel (all ten differed, canvas size included) — see §8.1 |
| `suppl/Supplementary_Tables_S1-S5.xlsx` | **no** | `openpyxl` writes `docProps/core.xml` with `dcterms:created` / `dcterms:modified`, and the ZIP entry headers carry the run time. Every cell is identical; only the container bytes differ. Two consecutive runs gave different MD5s. |

**Consequence for anyone verifying this deposit:** compare the **`.csv` tables**, not the
`.xlsx` workbook or a single archive checksum. A differing `.xlsx` hash does **not** indicate
a different result — the `.csv` files are authoritative and the `.xlsx` is a
reader-friendly repackaging of exactly those tables.

### 8.1 Correction — the PNG condition is the FreeType build, not a version pin (2026-10-21)

An earlier version of the table above said PNG bytes reproduce "if the Matplotlib version
matches". That is not sufficient, and a first diagnosis of this discrepancy was itself
wrong, so both the claim and the diagnosis are corrected here.

Re-running the figure scripts in a clean virtual environment built strictly from
`REQUIREMENTS.txt` produced ten PNGs that all differed from the deposit, **including their
pixel dimensions** (`fig2` 3360×1830 → 3366×1832). The Python versions were not the cause:
`matplotlib` was 3.10.6 in both environments and the bundled `DejaVuSans.ttf` was
byte-identical (`md5 49c0f03e…`). The difference is the **FreeType library Matplotlib links
against**:

| Environment | matplotlib | FreeType | Result |
|---|---|---|---|
| conda (the environment used for this deposit) | 3.10.6 | **2.13.3** | figures match the deposit byte-for-byte |
| PyPI wheel, fresh venv | 3.10.6 | **2.6.1** | all ten figures differ, canvas size included |

Layout under `bbox_inches="tight"` is computed from glyph metrics, so a different FreeType
build shifts the bounding box by a few pixels and the entire image changes. FreeType is
compiled into the Matplotlib wheel or supplied by the conda package — **no line in
`REQUIREMENTS.txt` can pin it**, which is why `pillow` turned out to be irrelevant
(it was pinned anyway, being the version that produced these files).

**Practical consequence:** treat figure hashes as environment-specific. A third party
should establish agreement on the **`.csv` / `.json` tables** (byte-exact, §8) and judge
figures by their content rather than their hash. Every number reported in the manuscript
reproduces exactly; only the rendered canvases are tied to the binary stack.

---

## 9. What a third party can and cannot re-run (measured 2026-10-21)

The deposit was re-run from a clean checkout in a fresh virtual environment, treated as if
the re-runner had no access to the developer machine.

| Layer | What was re-run | Outcome |
|---|---|---|
| Artefact | re-rendered `figs/` from the shipped JSON | **18/18 CSV byte-identical**; PNG bytes are environment-specific (§8.1) |
| Computation | each Stage re-run from shipped inputs | 15 of 20 scripts run clean once the cache in §9.1 is fetched; the other 5 need data listed in §9.2 |
| Data | starting from raw GEO downloads | the cohort cache is fetched by the helpers in §9.1; single-cell matrices are not covered |

### 9.1 Fetch the cohort cache before anything else

`cohort_hunt/geo_cache/` holds tens of MB of public GEO data and is **not** shipped. There
are two kinds of file in it and they are fetched differently.

**Series matrices** — one per cohort, fetched on demand by the parser itself, so a first run
just takes longer. To fetch the main cohorts ahead of time (the list is fixed inside the
script; GSE117606 is not in it and is fetched on demand instead):

```bash
python geo_dl.py
```

**Platform annotation** — not fetched automatically, and the reason most Stage B/C runs fail
on a fresh machine:

| Platform | Route | Note |
|---|---|---|
| GPL570, GPL4133 | `GPLxxx.annot.gz` | fetched automatically |
| GPL25373 (GSE117606) | `Homo_sapiens.gene_info.gz` from NCBI Gene | `.annot.gz` returns **HTTP 404** |

For GPL25373 the platform rows are Entrez IDs (`10000_at`), so `load_annot()` falls back to
`load_entrez_map()`, which needs the 5.2 MB NCBI file:

```bash
python geo_dl2.py URL:https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
```

Without it GSE117606 is dropped silently — `load_annot` returns `{}` and `load_entrez_map`
returns `{}` because the file is missing, so the cohort looks like it has no annotation. The
55.8 MB `GPL25373_family.soft.gz` also works (`python geo_dl2.py SOFT:GPL25373` then
`python gpl_soft_map.py GPL25373`) but is not needed for this cohort.

**Without the cache, most scripts fail loudly — `three_layer_check.py` used to be the
exception that failed silently.** Skipping a cohort it cannot load is correct, but the
script then wrote the canonical `three_layer_*.json` names anyway: a first run without the
cache replaced the complete archived result with a truncated copy (GSE4183 lost its
carcinoma group, so its AUC disappeared) while still exiting 0. Since v1.0.2 it writes
`three_layer_*.partial.json` and exits **2** whenever a cohort was skipped, so the archived
files are never overwritten. A non-zero exit with `.partial.json` output means the cache is
missing, not that the code is broken.

### 9.2 Stages that depend on data not in the repository

| Stage | Script(s) | Requires |
|---|---|---|
| B | `lcm_diag.py` | `cohort_hunt/geo_cache/GSE15960_series_matrix.txt.gz` |
| C | `three_layer_check.py`, `gse117606_paired.py`, `r3_sens_paneth.py`, `drift_diag.py` | the cache and the GPL25373 mapping from §9.1 |
| D | `sc161277_scan.py` | `sc161277/` (10x matrices for GSE161277) |
| E | `sc303070_compo*.py`, `sc303070_niche.py` | `sc303070/GSE303070_All_cells_metadata.csv.gz` |

Stage E cannot run from the repository alone. Its outputs
(`sc303070_composition_final.json`, `sc303070_niche.json`) are shipped, so every downstream
stage still runs.

With the §9.1 cache in place, Stage C reproduces the archived tables exactly: all three
cohorts load, `three_layer_pattern.json` is byte-identical, and `three_layer_results.json`
matches value for value including the GSE4183 carcinoma group.

### 9.3 Enrichr output is not byte-stable

`enrichr_dual.py` posts to the Enrichr web service, which mints a fresh `uid` per
submission. Every enrichment term, adjusted p-value and gene list reproduced exactly; only
`uid` differed. Compare all fields except that one.

### 9.4 Compare values, not bytes, for JSON

Two shipped JSON files can differ byte-for-byte between two correct runs:

| File | Why | What to compare |
|---|---|---|
| `enrichr_results.json` | new Enrichr job `uid` per submission | everything except `uid` |
| `three_layer_results.json` | the per-cohort `n` dict is built by iterating over a set of class labels, so key order varies with string hashing | parsed content |

Neither difference is a numerical one; equal JSON after `sort_keys=True` is the right test.
`three_layer_pattern.json`, `r3_sens_paneth.json` and `gse117606_paired_results.json` are
byte-identical across runs.

