## Promoter update: use the v2 module

See [R1M8_PROMOTER_README.md](R1M8_PROMOTER_README.md). The old audit_promoters.py and old archive retain historical results; use audit_promoters_v2.py for the latest promoter analysis. All three predefined comparisons remained non-significant. The 63-gene input has not been frozen for all subsequent M2/M6 work.

# Snapshot status and interpretation

## Version layers

1. `scripts/`, `results/` and `figures/` originate from the September 4 File S1
   package unless a later source is identified in `SOURCE_CODE_PROVENANCE.tsv`.
   Paths and line endings have been adapted for publication.
2. `build_key_clade_figures.py`, `plot_circular_phylogenies.py` and
   `export_phylogeny_itol_package.py` are the later plotting/export sources.
   `current_phylogeny/` contains the corresponding current Newick/display files.
   Figure 2 SVGs under `figures/` are historical outputs; do not substitute them
   for the current display trees or infer that every figure is publication-final.
3. `workflows/` and `supplement_20260908/` contain subsequent mapped/unassigned-read
   discovery, full short-read/Iso-Seq validation and family audits. These remain a
   separate supplement until final candidate adjudication and manuscript integration.

## Open scientific review items

- nsLTP classification rules are under re-audit. Historical scripts include
  member-specific/legacy-membership branches and an uncalibrated 50-bit similarity
  criterion. Cross-species application and the 8-cysteine motif interpretation
  require further review. Publishing the actual code preserves these decisions
  for inspection; it does not validate them or freeze the final catalogue.
- The mapped-read supplement examines 26,234,614 genome-mapped but
  transcriptome-unassigned fragments. The reported 62.81% uses all Salmon-unassigned
  fragments as its denominator, not all sequencing reads. Assembled ORFs and HMM
  hits remain candidates, not automatically new AMP genes.
- Full short-read and Iso-Seq workflow runs have recorded outputs. Candidate
  model/family decisions and the final manuscript update are still pending.
- Files named `*COMPLETE.PASS`, historical reports saying "final", and old decision
  summaries record the stage that generated them. They are not evidence that all
  reviewer requests or biological validation have been completed.

## Publication checks

Checksums and source-to-release provenance are provided. Publication preparation
checks syntax, the included synthetic workflow contract tests and archive integrity;
it does not rerun the full bioinformatics workflow. Consult `PUBLICATION_CHECKS.json`
for the actual checks and their results.
