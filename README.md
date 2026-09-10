**Final merged update (2026-09-10):** [Files S10/S11, current Tables S1–S59 and current figures](REVISION_UPDATE_20260910.md) include both the Tamarix anchor-input clarification and the [treatment-source clarification](REVISION_SUPPLEMENT_20260910.md). Use the merged workbook and merged checksums linked in the update guide, together with the preceding fixed release and 14-pattern promoter amendment.

> **Latest supplement (10 September 2026):** [Files S10–S11 and consolidated Tables S1–S59](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/tag/revision-20260910) add the characterization/duplication evidence and restored Arabidopsis comparison/current Figure 8. [Scope and treatment-source clarification](REVISION_SUPPLEMENT_20260910.md) distinguishes the conflicting published photoperiods; Table S25 has been corrected. Use this supplement together with the fixed 9 September release and the 14-pattern promoter amendment below.

> **Latest author figure formats (2026-09-09):** [Final Figure 3 and Figure 7](AUTHOR_FIGURE_UPDATE_20260909.md) provides the final author PDF/PNG/outlined-SVG layouts. This presentation-only update leaves the [14-pattern promoter amendment](PROMOTER14_UPDATE_20260909.md) and all numerical results unchanged.

> **Current analysis amendment (2026-09-09):** [The 14-pattern promoter and Figure 7 update](PROMOTER14_UPDATE_20260909.md) supersedes the earlier 15-pattern Figure 5/tables and supplies the current author Figure 7 with **Non-DE**. Retained counts and raw P values are unchanged; BH correction now uses 14/42 tests, with no significant enrichment. Read the fixed release below together with this amendment and the preceding reader-label update; archived releases are unchanged.

> **Current figure and reader-label amendment (2026-09-09):** See [the presentation update](READER_LABEL_UPDATE_20260909.md) for current A/R model-source labels, author figure layouts, explanatory table labels and portable plotting edits. The fixed scientific release below remains unchanged.

# R. soongarica AMP-family research: current revision

Current public revision: **2026-09-09**, tag **revision-20260909**.

[Download the versioned release](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/tag/revision-20260909) for the current scientific code,
processed evidence, source figures, Tables S1–S50 and recovered execution records.
The 2026-09-08 files and original `scripts/` / `workflows/` are historical sources;
their earlier catalogue status must not be read as the current result.

## Downloads

- [Complete current source and evidence bundle](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/download/revision-20260909/RsoAMP_current_sources_and_evidence_20260909.zip)
- [Current code only](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/download/revision-20260909/RsoAMP_current_code_20260909.zip) (data and model outputs are in the complete bundle)
- [Current source index](https://github.com/tarimteacher1/rsoamp-genes-revision/blob/revision-20260909/PUBLIC_SOURCE_INDEX_20260909.tsv)
- [Release checksums](https://github.com/tarimteacher1/rsoamp-genes-revision/blob/revision-20260909/RELEASE_SHA256_20260909.tsv)
- [Execution provenance](EXECUTION_PROVENANCE_20260909.json)
- [Run order and scope](RUN_ORDER.md)

## Current scientific scope

The primary R. soongarica catalogue contains 73 computational candidates:
15 defensins, 18 Snakin/GASA and 40 nsLTPs. RsLTP20 is a separate,
reference-dependent C4-like extended candidate. Seven reconstructed models are
distinct from six historical incomplete backfills. The current four-strategy
phylogenetic analysis reports six supported local groups, without a complete
formal subfamily classification. The 66 members represented in the original
RNA-seq reference yield 18 DE calls at S400 versus CK; seven reconstructed models
remain unquantified in that reference.

The T. austromongolica comparison uses 73/38 primary candidates, 27 reciprocal
best-hit pairs and 25 direct candidate pairs in the existing genomic-anchor
reference. The separate T. chinensis expression comparison uses a historical
63-query subset and an official-CDS-derived reference. Twenty pairs pass the
retained-hit/forward-coverage rules; 16 also have collinearity support. Each
contrast has 12 pairs with usable statistics in both studies. Four unique pairs
contribute seven concordantly upregulated, jointly significant comparison rows;
three pairs remain under the additional reverse-coverage rule. RsLTP17 is
coverage-sensitive. These results do not establish strict one-to-one orthology,
conserved function, an ecological-group trend or a between-species treatment effect.

## What was recovered in this release

The original successful v2 BLAST calls were recovered from retained execution
records: forward/reverse `-max_target_seqs 20/5`, `-evalue 1e-5`, `-seg yes`,
`-comp_based_stats 2`, eight threads and the standard 12 reported columns.
The recorded outputs contain 325/151 rows and match the archived result hashes.
Highest-score uniqueness remains restricted to retained hits. The MCScanX
invocation and successful completion record were also recovered. The executable
path is linked to Bioconda mcscanx 1.0.0 build h9948957_0 by the installed package,
binary hash and installation chronology; this is not a version printout or binary
digest collected during the original run. Recovery did not rerun the analyses.

## Package boundaries

Modules S2–S5 preserve the earlier outgroup, comparative, model-evidence and
localization resources. S6 contains the unified catalogue evidence. S7 contains
the independent T. chinensis expression reanalysis and annotation context. S8
contains the synchronized catalogue analyses and eight main figures. The updated
S9 contains the supported candidate-pair expression results and recovered
execution provenance. A small number of dated README/status statements in
historical modules describe their earlier stage; see STATUS.md for current mapping.

Scientific numeric results and sequence identifiers were not changed to prepare
the public distribution. Private execution prefixes were replaced with illustrative
paths. Manuscripts, reviewer correspondence, their document-production helpers,
compiled database indices and software installations are excluded. These
publication transformations and excluded files are indexed in the bundle. Source
tables and analytical/plotting code are retained. Public input accessions and
external dependencies are documented; this is not a preconfigured installation.
Use the publication manifest to verify the distributed bytes. PREPUBLICATION
manifests identify supplied originals before path normalization, not released bytes.

Independent qRT-PCR, peptide activity and causal salt-tolerance validation remain
absent. Computational support is not experimental validation. The existing
LICENSE_NOTICE.md continues to apply; no new software license is asserted.
