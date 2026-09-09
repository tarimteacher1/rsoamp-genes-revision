# Current run order and reproducibility scope

1. Download `RsoAMP_current_sources_and_evidence_20260909.zip` from [release revision-20260909](https://github.com/tarimteacher1/rsoamp-genes-revision/releases/tag/revision-20260909) and extract it into one directory.
2. Run `python verify_public_release.py` from that directory. It verifies the public
   file manifest using only the Python standard library.
3. To recount the existing paired-expression results, run
   `python modules/S9/verify_reported_counts.py`.
   This checks 16 eligible / 12 usable pairs per comparison, 3/4 jointly-DE pairs,
   and 6/9 same-sign effects at S200/S400. It is output verification, not a sequence search.
4. Consult the module README, execution records and source index before any
   scientific rerun. Configure `/path/to/...`, `${PAIRING_ROOT}`, `${EXPRESSION_ROOT}`,
   `${RSO_PROJECT}` and `${CONDA_ENVS}` for the local inputs and tools. Obtain the
   public assemblies, annotations, reads and third-party dependencies by their
   recorded source accessions; large raw reads, BAMs and installed databases are
   not redistributed. Existing paths do not imply a turnkey Windows workflow.
5. Use the earlier File S1 snapshot for the original genome/RNA-seq pipeline, then
   the versioned S6/S8 source code for catalogue reconciliation/current figures.
   Use S7 for T. chinensis RNA-seq and S9 for the existing candidate-pair comparison.
   The recovered S9 execution excerpts document the original order/settings and
   must not be relabelled as new analysis or unlimited-hit searches.

The code-only archive preserves relative module paths but does not include the
input data needed to execute analysis or plotting scripts. Use the complete bundle
for inspecting and verifying the retained results. No rerun of raw-read alignment,
phylogenetic inference or sequence searching was performed during publication.

Historical module checksums and PASS records describe their dated source files.
Public path normalization changes some byte hashes. The root
MANIFEST_SHA256_PUBLIC.tsv is authoritative for this release; S9's local manifest
was refreshed for its released files. It is expected that PREPUBLICATION checksums
do not match path-normalized bytes.
