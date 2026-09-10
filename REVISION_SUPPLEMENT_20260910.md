# Supplement to the 9 September public revision (10 September 2026)

Download the files from https://github.com/tarimteacher1/rsoamp-genes-revision/releases/tag/revision-20260910

- File S10: Tables S51–S56, precursor physicochemical properties, duplication classifications, seven within-species anchor pairs and nsLTP core-cysteine descriptors. Its current calculation script reproduces all six supplied TSVs from retained inputs. The original duplicate_gene_classifier command/version printout remains unavailable; upstream MCScanX was not rerun.
- File S11: Tables S57–S59, 26 Arabidopsis candidate-associated anchor pairs involving 22 red-sand members and 22 Arabidopsis loci; 13 red-sand members link to both comparator species. Source records, recovery/plotting scripts, the full neighborhood atlas and the current author Figure 8 are included. Original DIAMOND/MCScanX inference was not rerun. The table-recovery script reproduces the supplied TSVs.
- Supplementary_Tables_S1-S59_20260910.xlsx is the consolidated current workbook. The independent within-species check in S11 maps to formal Table S55 in S10, not a duplicate numbered table.
- S25_Treatment_source_clarification.tsv provides the current treatment provenance. No expression estimates or sample assignments changed during this correction.

## Treatment-source clarification

The original data paper (Song et al., 2024, Scientific Data 11:812, https://doi.org/10.1038/s41597-024-03644-y; Methods: Nucleic acid extraction and quality assessment) reports 26/16 degrees C, 16 h light/8 h dark, and half-strength Hoagland solution every three days. It identifies tender leaves, 0/200/400 mM Na2SO4 and three biological replicates per condition. Its three-week harvest explicitly concerns full-length transcriptome sequencing; it does not establish the age or treatment duration of the salt-treatment short-read samples.

Zhu et al. (2026, Chinese Bulletin of Botany 61:449–461, https://doi.org/10.11983/CBB25066; section 1.6, printed p.451) describes eight-week-old seedlings, one week of treatment and eight seedlings per replicate, and states that this follows the earlier transcriptome preparation. However, its reported photoperiod is 10 h light/14 h dark. These published accounts conflict. The actual photoperiod of the deposited RNA-seq samples has not been independently resolved; later-only details must not be labelled independently verified for those libraries. Table S25 distinguishes the sources and preserves the printed light-intensity unit as an unresolved source-unit issue.

## Earlier deposits and precedence

Use this supplement together with release revision-20260909 and PROMOTER14_UPDATE_20260909.md. The promoter analysis remains 14 patterns/42 tests with no significant enrichment. No previous release has been replaced. Historical workbook copies inside modules are retained solely as source inputs and do not supersede this S1–S59 workbook. File S11's author Figure 8 supersedes the earlier two-species display. The current manuscript's Ka/Ks discussion distinguishes 18 of 19 usable estimates below one from the 16 that also have nominal Fisher P < 0.05; the underlying table was already correct.

Private execution prefixes use illustrative paths. Manuscripts and reviewer correspondence are excluded. The released manifests identify distributed bytes; PREPUBLICATION manifests identify original retained records. Public distribution and TSV reproduction are not a new independent rerun of all biological analyses or experimental validation. The repository LICENSE_NOTICE.md continues to apply.
