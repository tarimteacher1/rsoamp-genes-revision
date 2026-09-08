# Expression reporting correction, 2026-09-04

This is a scoped update of the 2026-09-02 analysis archive. No DE model was refitted.
All 63 primary genes have TPM records. The original prefilter retains genes whose
total DESeq2 integer count across nine libraries is at least 10: 49 pass, 11 have
all-zero counts, and RsDEF6/RsLTP6/RsLTP21 have total counts of 3/2/6.
These 14 are annotated genes, not unquantified custom rescue models.

The revised reporting tables distinguish quantification, prefilter eligibility,
and contrast-specific adjusted-P availability. A prefiltered or adjusted-P-
unavailable result has an NA DE call, not FALSE. Six retained gene-contrast
records have a finite raw P value but no adjusted P value after independent
filtering; they are not part of the 14 prefiltered genes. The original numerical
estimates and the 7/17/12 positive DE sets are unchanged.

## Reproduce the audit

Use the original revision execution tree, with its Salmon quant.sf files,
rnaseq_sample_sheet.tsv, transcript_to_gene.tsv, frozen gene_TPM.csv, and
deseq2_dataset.rds. Raw/public inputs and the frozen model are not redistributed
by this archive. The audited environment is recorded in audit_R_sessionInfo.txt.
Use a new audit output directory; the auditor refuses to overwrite one.

```sh
Rscript scripts/audit_expression_prefilter.R /path/to/revision_R1_20260902 /path/to/new_prefilter_audit
python scripts/finalize_expression_catalogue.py   --catalogue /path/to/revision_R1_20260902/03_annotation_rescue/final_amp_primary_catalogue.tsv   --de-results /path/to/revision_R1_20260902/06_expression_reanalysis/results/all_gene_DE_results.csv   --tpm /path/to/revision_R1_20260902/06_expression_reanalysis/results/gene_TPM.csv   --normalized-counts /path/to/revision_R1_20260902/06_expression_reanalysis/results/gene_normalized_counts.csv   --vst /path/to/revision_R1_20260902/06_expression_reanalysis/results/gene_VST.csv   --prefilter-audit /path/to/new_prefilter_audit/primary_AMP_prefilter_audit.tsv   --output-dir /path/to/new_corrected_tables
python scripts/test_expression_status.py
```

The full run_expression_deseq2.R now exports the all-gene prefilter decision and
integer counts before applying the unchanged filter. This full analysis script
was not rerun for the reporting repair. Existing command logs are historical;
the reporting finalizer now requires the explicit --prefilter-audit argument.

## Version boundaries

The figures and other analyses retained in this archive are the frozen
2026-09-02 evidence assets. They are not replacements for the 2026-09-03
key-clade display figures embedded in the corrected manuscript. Historical
document-generation scripts are not a one-command build of the later manuscript.
Only the expression reporting tables, two expression scripts, the relevant old
document-generator wording, this audit, and the checksum manifest were updated.
Other review items, including AMP residual-read search and experimental
validation, are not resolved by this correction.
