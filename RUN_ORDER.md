# Workflow order and setup

## Before running

Unzip `RsoAMP_reproducibility_20260908.zip`. Its root contains `scripts/`, `results/`,
`current_phylogeny/`, `workflows/` and the supporting documentation.

Use Linux for the Bash and external bioinformatics workflows. Python scripts use
Python 3; dependencies vary by module (including Biopython, NumPy, pandas,
matplotlib and pysam). R modules use tximport and DESeq2. External tools include
STAR, Salmon, StringTie, gffread/gffcompare, HMMER, DIAMOND/BLAST, MAFFT, trimAl,
IQ-TREE, MCScanX, DeepSig, TMbed and PredGPI. Check the scripts you intend to run
and the recorded software inventory; not every dependency is needed for every step.

Replace `/path/to/rsoamp`, `/path/to/user-home` and any remaining illustrative
paths with your input/work directories and installed tool locations. Original
scripts often expect `revision_R1_20260902/<module>`; the packaged module outputs
are under `results/<module>`. Choose a new working copy and adapt `REVISION`,
`PROJECT`, command-line arguments and tool paths as each script supports them.
Do not run analysis directly against the retained source snapshot.

## Main workflow

| Order | Module | Main entry points |
| --- | --- | --- |
| 1 | Public input and read integrity | `download_verify_rnaseq.sh`, `download_verify_isoseq.sh`, `run_fastp_qc.sh` |
| 2 | Family curation | `build_nsltp_evidence.py`, `build_curated_nsltp_reference.py`, `finalize_nsltp_catalogue.py` |
| 3 | Four-strategy phylogeny | `prepare_phylogeny_inputs.py`, `run_phylogeny_sensitivity.sh`, `summarize_phylogeny_sensitivity.R` |
| 4 | Historical backfill and catalogue | `audit_backfill_loci.py`, `evaluate_isoseq_gene_models.py`, `audit_genome_wide_rescue.py`, `build_final_amp_catalogue.py` |
| 5 | Genome mapping and unassigned reads | `align_rnaseq_genome.sh`, `run_unmapped_read_classification.sh`, `run_unmapped_amp_search.sh` |
| 6 | Comparative genomics | `build_tamarix_amp_evidence.py`, `run_tamarix_nsltp_recuration.sh`, `run_rso_tamarix_mcscanx.sh`, `build_comparative_outputs.py` |
| 7 | Expression | `run_salmon_requant.sh`, `run_expression_deseq2.R`, `audit_expression_prefilter.R`, `finalize_expression_catalogue.py` |
| 8 | Promoter and Ka/Ks audits | `audit_promoters.py`, `audit_kaks.py` |
| 9 | Figures and tables | `plot_revision_*.py`, `build_supplementary_workbook.mjs`; latest tree displays: `build_key_clade_figures.py` |

Historical result/command records describe what was run. This table is an index,
not a substitute for input preparation or evidence that the family rules are settled.
The expression correction is detailed in `EXPRESSION_STATUS_CORRECTION_20260904.md`.

## Later mapped-read supplement

`workflows/mapped_rescue/Snakefile` calls `run_stage.py` for preparation, extraction,
STAR mapping, assembly, screening and evidence aggregation. Adapt
`config.example.json`, including its `configuration_file` field. Run Snakemake
from a new run directory using that absolute config and Snakefile path.
The source/config directory and the output run directory should be distinct.

Then follow the dependent modules `full_read_validation`, `overlap_family_audit`
and `recursive_family_audit`. Their shell scripts declare input and output paths
near the top. These workflows require large external inputs absent from this
repository. Later candidate evidence is not integrated into the September 4 catalogue.

## Small verification without research data

From the extracted package root:

```bash
python workflows/mapped_rescue/test_contract.py
python verify_manifest.py
```

The first command uses synthetic data to check read pairing, interval and ORF
logic. The second verifies the files shipped in this snapshot. Neither command
tests the full biological workflow or independently validates the scientific claims.
