# Genes RsoAMP R1.M6 focused supplement v1.0.0

Native Snakemake wrapper for the project's STAR/StringTie/AMP evidence workflow.
No input downloads, installations, background services, deletion, DE refitting,
automatic primary-catalogue changes, or publication-ready claim generation.

## Input contract

An absolute JSON config contains `revision`, `genome`, `hmm_models`, `tools`
(absolute paths for python, STAR, samtools, stringtie, gffread, gffcompare,
diamond, hmmsearch, getorf), and `samples` (run, treatment, replicate).
These point to the verified remote original execution tree. Inputs and config
are outside this immutable source directory. Outputs are relative to `results`
inside a new run directory, never inside the frozen revision.
The `configuration_file` field repeats that same absolute JSON path, for the
stage runner; it is bound by the native plan rather than guessed from cwd.

## Native entrypoint

`Snakefile` calls `run_stage.py` with stage-specific actions and the immutable
config file. Rules use a total resource cap and preserve exit codes. The runner
uses explicit installed binary paths and records every external command.
Each tool's bin directory is prepended to PATH for that child process only,
including the existing EMBOSS wrapper's `_getorf` and data-directory lookup.

## Stages

prepare -> extract paired target fragments -> per-sample STAR -> per-sample
StringTie -> merge/gffcompare/gffread -> getorf/HMMER/DIAMOND -> evidence audit.
All 9 sample alignments/assemblies are retained. A successful discovery marker
does not mean complete local transcript validation or closure of reviewer M6.

## Scientific limits

The target is the 26,234,614 genome-mapped but transcriptome-unassigned fragments,
not the 5,863,532 decoy-assigned fragments or all raw RNA-seq. Protein references
cover the three retained AMP families; the HMM scan covers six historical
families. Stop-delimited ORFs >=30 aa and assembled transcripts >=90 nt define
an explicit sensitivity floor. Broad prolamin similarity is not an nsLTP call.
Candidates remain provisional pending full-read local evidence, secretion and
family-boundary validation. No activity or novel gene claim is made.

`unique_fragments` counts distinct fragment IDs with at least one overlapping,
strand-compatible primary alignment whose NH is 1; `multimapping_fragments`
counts IDs with an overlapping primary NH>1 alignment. These can overlap if
mates differ. Coverage uses only NH=1 aligned blocks, never intron spans.
Known-AMP gene-body/flank coverage is region coverage, not coding-exon coverage.

## License

Project analysis code; no public redistribution license was assigned by the
authors. Installed third-party tools and databases are not redistributed.

## Verification

`python test_contract.py` checks read naming, paired FASTQ validation, count
reconciliation, GTF interval logic and ORF boundary parsing with synthetic data.
Read-only readiness checks do not execute the scientific workflow.
