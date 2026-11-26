# Release notes — v1.0.5

Reproducibility package accompanying the manuscript submitted to *BMC Genomics*:

> **Cellular composition and epithelial cell state contribute in opposite directions to the
> colorectal adenoma–carcinoma transcriptome**

## Changes in v1.0.5 (minor-revision alignment; robustness analyses deposited)

This release aligns the archive with the revised manuscript (minor revision round) and, for
the first time, ships the code behind every robustness and sensitivity number that the
revised Methods section reports.

**New in the deposit — Stage J (robustness and sensitivity):**

- `_supp_bootstrap_auc_ci.py` — bootstrap CIs (B = 3000) for the headline AUCs.
- `_supp_permutation.py` — sample-level label permutation of the component means (the test
  the component-level conclusions now rest on).
- `_cam_gene_correlation_test.py` — module-block sign-flip and CAMERA-style variance
  inflation, reported alongside the sample-level permutation. Reads the module assignment
  newly shipped under `upstream_r/02.WGCNA/gene_module_assignment.csv`.
- `_supp_doublet.py` — doublet sensitivity of the single-cell decomposition.
- `_supp_sva_sensitivity.py` — SVA-lite orthogonalisation sensitivity, with its figure.
- `_loo_macro_sensitivity.py` — leave-one-cohort-out stability of the consensus signature.
- `_gse8671_validation.py` and `_dl_gse8671_v5_soft.py` — independent validation on
  GSE8671 (32 paired endoscopic biopsies). The series-matrix copy served by NCBI is
  corrupted at a fixed byte offset; the helper reads the SOFT family file instead.
- `_m1_drift_arithmetic_sim.py` — arithmetic simulation showing that the sample-level
  drift correction absorbs any uniform offset, with its figure.
- Archived outputs of all of the above ship as JSON; the two robustness figures sit at the
  repository root (they are not manuscript additional files).
- `GSE8671` and `GSE109308` are added to the public-input table (§2 of
  `REPRODUCIBILITY.md`), which previously covered ten of the twelve series.

**Supplementary-table numbering re-mapped.** The revised manuscript re-labelled its
supplementary tables so that they follow first-citation order. The shipped file names and
worksheet names now carry the manuscript's numbers:

| manuscript table | shipped files |
|---|---|
| Table S3 (a/b) | `S3a_sc161277_sample_QC.csv`, `S3b_sc161277_decomposition.csv` |
| Table S4 | `S4_four_cohort_performance.csv` |
| Table S5 (a/b) | `S5a_panel_pattern_matrix.csv`, `S5b_paneth_sensitivity.csv` |

Cell content is unchanged — every re-named file is byte-identical to its previously shipped
counterpart; only the names moved.

**Manuscript-build scripts updated.** `build_bmc_manuscript.py` follows the revised
abstract, the new "Robustness and sensitivity" Methods section, and the re-mapped
additional-file list; `_bmc_references.json` carries 32 references (one added: the
GSE8671 source publication). `make_figs.py` refinements: Figure 2 axis label and
annotations now state the pair-level proportions directly (53/56 and 8/56), and Figure 3's
component ratio is written as 0.316×. `REPRODUCIBILITY.md` gains a Stage J section and a
12-series input table.

No reported number changed value in this round except as documented in the manuscript
revision itself; the archived analysis JSONs for Stages A–I are byte-identical to v1.0.4.

## Changes in v1.0.4 and earlier

See `RELEASE_NOTES_v1.0.4.md`, shipped alongside this file.
