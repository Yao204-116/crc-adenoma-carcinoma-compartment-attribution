# Release notes — v1.0.4

Reproducibility package accompanying the manuscript submitted to *BMC Genomics*:

> **Cellular composition and epithelial cell state contribute in opposite directions to the
> colorectal adenoma–carcinoma transcriptome**

## Changes in v1.0.4 (title, alignment, one figure)

**This release exists because a published title cannot be edited.** The manuscript title
changed after v1.0.3 was archived, and a Zenodo record freezes its metadata at publication.
Anyone following the DOI printed in the *Code availability* statement therefore arrived at a
record describing a differently titled study. The concept DOI (10.5281/zenodo.23498178) always
resolves to the newest version, so the fix is a new version under the same concept, not a new
record.

While preparing it, three accumulated differences between the shipped code and the code that
produced the submitted manuscript were closed:

- **Every GEO series now carries its source publication.** `_bmc_references.json` grew from 21
  to 31 entries. The association was read from the `Series_pubmed_id` field of each series'
  SOFT record rather than from the PubMed summary, which returns an empty list for most
  series and had previously been taken to mean that no paper was linked. Seven series have an
  official link. For the remaining four the citation rests on the submitter and the paper
  authors being the same group, on a matching platform and tissue type, and on the accession
  appearing in the paper itself; those four are marked as inferred in the internal provenance
  table, not presented as official.
- **`verify_bmc.py` grew from 24 to 31 checks.** One of the new ones fails if any accession
  first appears without a citation bracket. It matches `GSE\d{4,6}`, not `\d{5,6}`, because
  `GSE4183` has four digits and a five-or-more pattern skips it silently — the check would
  have reported a clean pass over ten series instead of eleven.
- **The panel heat map is no longer folded.** `make_figs.py` used to plot
  `max(AUC, 1 − AUC)`, which turned a component separating the stages in the opposite
  direction (AUC 0.111) into 0.889. The legend asked the reader to look for values below 0.5,
  and no such value could appear. Cells are now plotted as computed on a 0–1 scale, so the
  reversal the Results describe is visible in the figure.

Scoring is now stated the same way everywhere: the undirected mean over signature genes
(defined in the Methods) is the primary number, and the directed difference between the two
component means appears alongside it in the text. `qa_numbers.py` checks the undirected values.
`build_bmc_manuscript.py` no longer numbers the section headings, matching the BMC style.

**No reported number, table or figure other than the heat map differs from v1.0.3.** Nothing
in the analysis pipeline changed; the data, the models and the statistics are the same.

The deposit is timed to the revised manuscript rather than the first submission: the author
list now lists Chengkun Yao and Bo Jiang as co-first authors and Wanting Li as third author,
the abstract follows the four-part Background/Methods/Results/Conclusions structure, and the
submitted text was shortened and re-edited. None of these edits touch the analysis, so the
scripts and derived tables here are unchanged in content apart from the heat-map fix above.

## Changes in v1.0.3 and earlier

See `RELEASE_NOTES_v1.0.3.md`, which is shipped alongside this file. v1.0.3 corrected one
README line (a doubled `https://doi.org/` prefix). v1.0.2 was a reproducibility fix following
an external audit: developer-local absolute paths removed, missing inputs shipped, and the
sample classifier no longer degrading silently to substring matching.
