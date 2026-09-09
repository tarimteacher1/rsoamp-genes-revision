# 14-pattern promoter and Figure 7 update — 2026-09-09

Download [the amendment ZIP](RsoAMP_promoter14_figure7_20260909.zip) and verify its
[SHA-256](RsoAMP_promoter14_figure7_20260909.zip.sha256).

This update supersedes the 15-pattern promoter analysis in the fixed revision-20260909 release
and reader-label package. It contains the current 14-column Figure 5, source code, cached
counts/cohort records, and replacement Tables S9/S10/S27/S29. Historical 15-pattern records are
retained explicitly for auditing; old releases remain unchanged. Existing reader-label material
continues to apply to all other figures and tables. The included S28 is only an unchanged historical
reference and must not overwrite the corrected current S28 reader labels.

The unsupported custom BACGTG pattern was removed on provenance and interpretability grounds,
not statistical significance. BACGTG includes G-box when B=C and contains the ABRE-associated
ACGTG core. Remaining motifs can still overlap; no functional independence is implied.

All 305,172 retained gene-pattern counts and 42 raw one-sided Fisher P values are unchanged.
Within-scenario BH now uses 14 tests; joint BH uses 42. Minimum BH14 values are 0.913467,
0.593905 and 0.719358 for the primary, clean, and clean-with-UTR scenarios. All joint BH42
values are 1.0. No motif is significantly enriched at BH FDR < 0.05. Existing sequences and
boundaries were retained; no new sequencing, genome extraction, or expression analysis was done.

The current Figure 7 author PDF/PNG uses **Non-DE** for grey points: genes not meeting the
joint BH-adjusted P < 0.05 and absolute unshrunken log2 fold change > 1 criteria. Missing or
unquantified members are not recoded as Non-DE. Its numerical data and 8/18 DE calls are unchanged.
The older reader-package SVG with a long grey-point label is not the final author Figure 7 display.

Inside the ZIP, run `python verify_patch.py`, then `cd promoter14` and `python run_all14.py`.
The package supplies exact count, raw-P and adjusted-P comparison records. It contains research
data, code, figures and their legends, and excludes manuscript files and reviewer correspondence.
