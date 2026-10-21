# Release notes — v1.0.2

Reproducibility package accompanying the manuscript submitted to *BMC Genomics*:

> **Composition or cell state? Decomposing the adenoma–carcinoma transcriptome across
> tissue, epithelium and single cells**

## Changes in this release

**v1.0.2 is a reproducibility fix. No number reported in the manuscript, no table and no
figure changed.** It follows an audit in which the release archive was downloaded and re-run
from a clean environment, as an outside reader would. One shipped intermediate
(`drift_diag.json`) is refreshed so that it matches what the code now produces; the group
that appears in it is read by nothing (§ below).

What the audit found and what changed:

- **Developer-local paths removed.** `sig_validate.py` and `gse117606_paired.py` resolved
  their working directory from a hard-coded absolute path on the development machine,
  contradicting the README's statement that paths are resolved relative to the script.
  Both now use `os.path.dirname(os.path.abspath(__file__))`. On a third-party machine the
  old code silently read nothing and produced truncated results with exit code 0.
- **Inputs that stages read but that had never been shipped are now included.**
  `consensus_signature_tier.csv`, `manuscript_facts.py`, `epithelial_core_genes.json`,
  `stromal_layer_genes.json`, `cohort_hunt2.py`, `gpl_soft_map.py` and
  `collide_results.json`. `enrichr_dual.py` crashed at import time without the two JSON
  files; the enrichment results it now reproduces match the archived `enrichr_results.json`
  in every term and p-value (the Enrichr job `uid` differs by design). `cohort_hunt2.py` is
  the omission that would **not** have announced itself: `sig_validate.py` wraps
  `from cohort_hunt2 import classify` in `try/except` and falls back to naive substring
  matching. Measured on the v1.0.1 archive, that fallback classified GSE4183's 15 carcinoma
  samples as `other`, so the cohort's AUC disappeared and `three_layer_results.json` was
  written truncated — with exit code 0. A wrong answer that runs is worse than an
  `ImportError`.
- **`sig_validate.py` no longer erases the four-cohort results when run without
  arguments.** Its main loop was `for gse in sys.argv[1:]`, so `python sig_validate.py`
  iterated zero cohorts and then dumped `{}` over `cohort_hunt/sig_validation_results.json`
  — the file `gse117606_paired.py` backs up and `manuscript_facts.py` turns into
  `FACTS['four_cohort_auc']`, which `build_suppl_tables.py` (the README's smoke test) reads.
  It now defaults to the four cohorts it is named after and refuses to write at all when
  nothing was computed. The main loop also only runs as `__main__`, because three stages
  load this file with `exec_module` purely to borrow its functions.
- **`collide_results.json` shipped, and the four-cohort fallback widened.** That file is
  produced by `collide_all.py` and read by `manuscript_facts.py`, but was never packaged —
  so re-running `manuscript_facts.py` yielded a version without the `collide` block.
  `manuscript_facts.py` additionally now prefers `cohort_hunt/sig_validation_results.json`
  and only falls back to the `_4cohort` backup, instead of depending on which script ran
  first.
- **`three_layer_check.py` no longer overwrites complete results with partial ones.** The
  cohort loader silently skips any cohort whose local `geo_cache/` is missing, but the
  script still wrote the canonical output filenames — so a reader's first run replaced the
  archived result with a truncated copy (GSE4183 lost its carcinoma group, so its AUC
  vanished, while the exit code stayed 0). It now writes `*.partial.json` and exits 2
  whenever a cohort was skipped, leaving the archived files intact.
- **`drift_diag.py` exec namespace fixed.** After the path change above, this script raised
  `NameError: __file__` on every machine, because it `exec`s the first 266 lines of
  `sig_validate.py`, which resolves its own directory from `__file__`. The traceback pointed
  at `<string>`, giving no hint that the code came from `exec`.
- **`drift_diag.json` regenerated under the shipped code.** The v1.0.1 copy was written while
  the classifier was degraded: its `GSE4183` block has only `normal` and `adenoma`, while a
  run of the shipped scripts also yields `carcinoma`. Only `GSE15960` is read by any figure
  (`make_figS1.py`), and that block is unchanged; no reported number moves. The file is
  replaced so that the archive matches what the archive's own code produces, rather than
  documenting a permanent discrepancy.
- **`REQUIREMENTS.txt` pins Pillow and corrects the reproducibility claim.** The earlier
  text said PNG bytes reproduce "if the Matplotlib version matches". Measured on a clean
  install, that is not sufficient: figure layout depends on the FreeType build that
  Matplotlib links against, which no pip pin can control. §8.1 of `REPRODUCIBILITY.md`
  records the measurement; §8 now tells readers to compare CSV/JSON for exact agreement
  and to treat figure hashes as environment-specific.
- **Python floor corrected, 3.9+ → 3.11+.** The pinned `numpy==2.3.5` does not install on
  3.9 or 3.10, so the previous statement could not be followed as written.
- **`REPRODUCIBILITY.md` §9 rewritten from measurement.** §9 now separates the two kinds of
  cache file and gives the working route for the one platform whose annotation is not on
  GEO (`GPL25373` → NCBI `Homo_sapiens.gene_info.gz`), states which stages need which data,
  and records that two JSON files can differ byte-wise without differing numerically
  (Enrichr `uid`; the `n` dict's key order in `three_layer_results.json`).

## What is in this release

| Component | Contents |
|---|---|
| `*.py` | Stages A–I of the analysis, the data-download helpers, and the shared modules the stages import (per-file inventory in `_MANIFEST.txt`) |
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

Python 3.11+ (developed on 3.13.9), see `REQUIREMENTS.txt`. R (limma, WGCNA) was used
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
