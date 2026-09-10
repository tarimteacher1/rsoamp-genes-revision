# Revision update: S10/S11, current tables and figures (2026-09-10)

Download the [10 September update release](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/tag/revision-20260910). It complements, rather than modifies, the preceding fixed scientific release and the 14-pattern promoter amendment.

- Full new source/evidence bundle: `RsoAMP_revision_additions_20260910.zip` in the release assets.
- [New analytical code](RsoAMP_added_code_20260910.zip); data are in the full bundle.
- [Current Supplementary Tables S1–S59](Supplementary_Tables_S1-S59_20260910.xlsx).
- [Source-to-public mapping](SOURCE_TO_PUBLIC_20260910.tsv), [full bundle content checksums](ADDITION_CONTENTS_SHA256_20260910.tsv), and [download checksums](RELEASE_SHA256_20260910.tsv).

Files S10 and S11 correspond to Tables S51–S56 and S57–S59 respectively. S10 covers physicochemical properties, duplication modes, seven within-species candidate anchor pairs and nsLTP structural descriptors. S11 retains the Arabidopsis comparison (26 primary candidate anchors; 22 Rso members and 22 Arabidopsis loci), integrates it with Tamarix in Figure 8, and provides all 26 Arabidopsis neighborhoods. Thirteen Rso members have links to both comparators. Current author PDFs and captions accompany the full bundle.

The current S21 scope wording makes the actual Tamarix anchor input explicit: 21,791 Rso reference records + 22,374 Tamarix reference records + five homology-derived Tamarix defensin models = 44,170. The seven new Rso models were not in this anchor input. The separate 44,177-protein sequence database includes those seven models. Scientific counts and results are not changed by this clarification or by publication.

The [9 September release](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/tag/revision-20260909) supplies the preceding S2–S9 evidence; the earlier repository snapshot supplies S1. Read it with the [reader-label update](READER_LABEL_UPDATE_20260909.md), [14-pattern promoter amendment](PROMOTER14_UPDATE_20260909.md) and this update. The authoritative current workbook is S1–S59 above; workbooks inside source modules are retained analysis inputs.

Full-bundle README gives dependencies, run order and verification. The S10 script recalculates added properties and checks retained duplication/anchor evidence. The S11 scripts re-parse retained DIAMOND/MCScanX outputs and redraw figures; they do not rerun those upstream tools. Source scripts, original logs and their provenance are retained where available. Missing historical version/command printouts remain explicitly identified; publication does not fill them with invented settings.

Private execution prefixes are normalized and mapped. No manuscript, reviewer correspondence, credentials, raw sequencing reads or installed software are distributed. Existing licensing notices remain applicable. Computational family and expression evidence does not establish antimicrobial activity or causal salt tolerance; no independent qRT-PCR validation is added.
