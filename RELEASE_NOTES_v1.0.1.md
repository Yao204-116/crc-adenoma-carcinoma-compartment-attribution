# Release notes — v1.0.1

Reproducibility package accompanying the manuscript submitted to *BMC Genomics*:

> **Composition or cell state? Decomposing the adenoma–carcinoma transcriptome across
> tissue, epithelium and single cells**

## Changes in this release

- Manuscript title updated to the version submitted to *BMC Genomics*.
- `GSE109308` (familial adenomatous polyposis / sporadic carcinoma single-cell series)
  added to the documented datasets: it is used in the Discussion as an independent
  boundary-condition check, not as an input to any effect-size estimate.
- Author metadata and funding statement added.

The analysis code and all derived results are unchanged from v1.0.0.

## What is in this release

| Component | Contents |
|---|---|
| `*.py` (34 scripts) | Stages A–I of the analysis, plus the data-download helpers |
| `*.json`, `*.csv` | The derived intermediate results — the numbers behind every figure and table |
| `02.result/` | Per-cohort DEG tables for GSE37364, GSE20916, GSE41657 |
| `upstream_r/` | The five CSV products of the upstream R pipeline (limma, WGCNA, enrichment) |
| `figs/` | Rendered main figures (1–8) and supplementary figures (S1–S3) |
| `suppl/` | Supplementary Tables S1–S5 (`.csv` and `.xlsx`) |
| `README.md` | Environment, layout, quick-start |
| `REPRODUCIBILITY.md` | Full script → intermediate → manuscript-artefact mapping |

## Data

No new primary data were generated. All input data are public (NCBI GEO):

| Accession | Role |
|---|---|
| GSE41258, GSE37364, GSE20916, GSE41657 | Discovery cohorts (whole tissue) |
| GSE4183, GSE117606 | Three-layer validation (whole tissue) |
| GSE15960 | LCM pure-epithelium cohort |
| GSE161277 | Single-cell cohort (normal → adenoma → carcinoma) |
| GSE109308 | Single-cell cohort (FAP / sporadic carcinoma) — Discussion boundary condition only |
| GSE128435 | Screened and excluded (composition-confounded) — retained for transparency |
| GSE303070 | CosMx spatial atlas — used only as a descriptive counter-example |

## Environment

Python 3.9+ (developed on 3.13), see `REQUIREMENTS.txt`. R (limma, WGCNA) was used
upstream; its outputs are shipped under `upstream_r/`.

## Reproduce

```bash
git clone https://github.com/Yao204-116/crc-adenoma-carcinoma-compartment-attribution
cd <repo>
pip install -r REQUIREMENTS.txt
python build_suppl_tables.py     # smoke test: rebuilds Supplementary Tables S1-S5
```

## License

MIT — see `LICENSE`.

## Known limitations

See `README.md` §6 and `REPRODUCIBILITY.md`: the GSE128435 exclusion and the
"non-arbitrable reference" status of GSE303070 are documented judgement calls rather than
automatic filters, and all p-values are reported as agreement between two exchangeable
calibers rather than as single significance tests.
