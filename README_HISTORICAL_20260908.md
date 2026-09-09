> **Promoter update (2026-09-08):** The latest promoter analysis is in [R1M8_PROMOTER_README.md](R1M8_PROMOTER_README.md) and [R1M8_promoter_revision_20260908.zip](R1M8_promoter_revision_20260908.zip). It supersedes the historical promoter module below: 15 double-stranded patterns, full-genome interval rescan, symmetric QC and UTR sensitivity checks. No enrichment in any predefined scenario. Other modules retain the version boundaries described below.

# R. soongarica AMP analysis code — Genes revision

Research code and selected processed outputs for the identification, curation,
phylogenetic analysis, comparative genomics and expression analysis of
*Reaumuria soongarica* AMP-family candidates.

**Snapshot: 2026-09-08. This is a revision-stage research snapshot.**
Read [STATUS.md](STATUS.md) before interpreting the catalogues or candidate calls.

## Start here

- [RUN_ORDER.md](RUN_ORDER.md): workflow order, path adaptation and a small verification command.
- [scripts/](scripts/): analysis and plotting scripts underlying the current working manuscript.
- [workflows/](workflows/): later mapped-read discovery and full-read/family follow-up workflows.
- [DATA_SOURCES.md](DATA_SOURCES.md): public input accessions and dependency records.
- `RsoAMP_reproducibility_20260908.zip`: downloadable scripts, selected processed results,
  decision tables, original ML trees, current display trees and source-version records.
- `MANIFEST_SHA256.tsv`: checksums for the files in this repository.

The working catalogue in the September 4 package contains 63 annotated genes:
9 defensins, 18 Snakin/GASA and 36 canonical nsLTPs (27 core and 9 predicted LTPg).
All 63 have TPM records; 49 pass the expression prefilter. The six historical
genome-backfilled ORFs are outside this primary catalogue.

The release adds the subsequent mapped/unassigned-read and full-read validation
workflow sources. Those candidate outputs have not yet been incorporated into a
final frozen catalogue. Computational evidence does not establish peptide activity;
independent qRT-PCR and peptide activity assays are not supplied.

## Reproduction scope

This repository documents the research execution history. It is not a turnkey
installation or a claim that every analysis was rerun during publication preparation.
Large public reads, genome assemblies, BAM files, indices, installed software,
third-party databases, manuscripts and reviewer correspondence are excluded.
Historical absolute project and home paths were replaced with `/path/to/...`
placeholders, and text line endings were normalized. Scientific algorithms and
thresholds were not changed for this publication. Configure the paths and obtain
the external inputs before a scientific rerun.

See [LICENSE_NOTICE.md](LICENSE_NOTICE.md) for the current permission status.
For citation in the manuscript or response letter, use this repository's actual
URL together with the specific Git commit and the snapshot date. No DOI is assigned
by this package.
