# Data sources for the current 2026-09-09 revision

Current module records and release provenance take precedence over dated historical
scope statements. Original-source accessions and checksums are retained with each
module in the public bundle.

- R. soongarica genome/annotation: Figshare 10.6084/m9.figshare.25533064.v2 and
  GenBank JBEBFM000000000; RNA-seq/Iso-Seq PRJNA1063761, SRR27540875–SRR27540884.
- T. austromongolica genome: GCA_039764185.1 and Figshare 10.6084/m9.figshare.25106726.
- T. chinensis genome/annotation: GigaDB 10.5524/102417. RNA-seq PRJNA855335;
  six selected runs SRR19973985/986/987/991/993/994 are documented individually in S7.
- Source gene/transcript identifiers, official-CDS reference checks, selected
  isoforms, software versions, original-call excerpts and reporting thresholds are
  preserved in S6–S9 and EXECUTION_PROVENANCE_20260909.json.

Private execution prefixes are illustrative placeholders. Full raw reads, software
installations and compiled third-party databases are not bundled. Obtain them from
the documented sources; retain their provider terms and version identifiers.

## Historical source note

# Public data and dependency sources

- *R. soongarica* genome paper: https://doi.org/10.1038/s41597-024-03644-y
- Genome/annotation deposit: https://doi.org/10.6084/m9.figshare.25533064.v2
- RNA-seq/Iso-Seq project: NCBI BioProject **PRJNA1063761**, SRA **SRP483612**.
  Exact run identities and checksums are in `results/00_audit/`.
- *T. austromongolica* genome: https://doi.org/10.1093/dnares/dsae021
  and https://doi.org/10.6084/m9.figshare.25106726
- SILVA 138.2: https://www.arb-silva.de/documentation/release-1382/
- *R. soongarica* chloroplast: **NC_041273**.
- *Myricaria laxiflora* mitochondrial proxy: **MW971331**.

These identifiers are transcribed from the execution package. Original input
checksums, run metadata, reference manifests and software versions are retained
in the downloadable archive, particularly `results/00_audit/software_inventory.tsv`.
The supplement contains its own effective configuration and software/input identity
records. No raw-read or third-party database redistribution is implied.

Historical per-module checksum records refer to the original execution inputs.
Because publication copies have normalized paths and line endings, use the root
`MANIFEST_SHA256.tsv` to verify this release. Do not treat an original absolute-path
checksum list as the checksum manifest of the edited publication copies.
