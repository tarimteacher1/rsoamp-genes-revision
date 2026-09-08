#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
revision <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
if (dir.exists(out)) stop("Audit output already exists; use a new directory")
dir.create(out, recursive=TRUE)
expr <- file.path(revision, "06_expression_reanalysis")
samples <- read.delim(file.path(expr, "rnaseq_sample_sheet.tsv"))
samples$condition <- factor(samples$condition, levels=c("CK", "S200", "S400"))
rownames(samples) <- samples$run
files <- file.path(expr, "salmon_quant", samples$run, "quant.sf")
names(files) <- samples$run
stopifnot(all(file.exists(files)))
tx2gene <- read.delim(file.path(expr, "transcript_to_gene.tsv"))
stopifnot(!anyDuplicated(tx2gene$TXNAME))
txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="no")
dds0 <- DESeqDataSetFromTximport(txi, colData=samples[,c("condition","replicate")], design=~condition)
count_mat <- counts(dds0)
keep <- rowSums(count_mat) >= 10
old <- readRDS(file.path(expr, "results/deseq2_dataset.rds"))
stopifnot(identical(rownames(count_mat)[keep], rownames(old)))
stopifnot(identical(colnames(count_mat), colnames(old)))
stopifnot(identical(unname(count_mat[keep,,drop=FALSE]), unname(counts(old))))
old_tpm <- read.csv(file.path(expr, "results/gene_TPM.csv"), check.names=FALSE)
old_tpm_mat <- as.matrix(old_tpm[match(rownames(txi$abundance), old_tpm$gene_id), samples$run])
stopifnot(max(abs(old_tpm_mat - txi$abundance)) < 1e-7)
status <- ifelse(keep, "RETAINED_TOTAL_COUNT_GE_10",
                 ifelse(rowSums(count_mat) == 0, "FILTERED_ALL_ZERO_COUNTS", "FILTERED_TOTAL_COUNT_LT_10"))
audit <- data.frame(
  gene_id=rownames(count_mat),
  tximport_estimated_count_sum=rowSums(txi$counts),
  deseq2_rounded_count_sum=rowSums(count_mat),
  nonzero_count_samples=rowSums(count_mat > 0),
  nonzero_tpm_samples=rowSums(txi$abundance > 0),
  max_TPM=apply(txi$abundance, 1, max),
  prefilter_threshold_total_count=10L,
  retained_for_DE=keep,
  prefilter_status=status,
  check.names=FALSE
)
write.table(audit, file.path(out, "all_gene_prefilter_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
catalogue <- read.delim(file.path(revision, "03_annotation_rescue/final_amp_primary_catalogue.tsv"))
idx <- match(catalogue$gene_id, audit$gene_id)
stopifnot(!anyNA(idx))
primary <- cbind(catalogue[,c("member_id","family","gene_id")], audit[idx,setdiff(names(audit),"gene_id")])
write.table(primary, file.path(out, "primary_AMP_prefilter_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene_id=rownames(count_mat), count_mat, check.names=FALSE),
            file.path(out, "all_gene_prefilter_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(member_id=catalogue$member_id, gene_id=catalogue$gene_id,
                       count_mat[match(catalogue$gene_id,rownames(count_mat)),], check.names=FALSE),
            file.path(out, "primary_AMP_prefilter_counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
source_files <- c(files, file.path(expr, "transcript_to_gene.tsv"), file.path(expr,"rnaseq_sample_sheet.tsv"),
                  file.path(expr,"results/deseq2_dataset.rds"), file.path(expr,"results/gene_TPM.csv"))
write.table(data.frame(file=source_files,md5=unname(tools::md5sum(source_files))),
            file.path(out,"source_md5.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"audit_R_sessionInfo.txt"))
writeLines(c("Retained gene IDs and integer counts match the frozen DESeq2 dataset exactly.",
             "Reimported gene TPM matches the frozen TPM matrix within absolute tolerance 1e-7.",
             "No DESeq model was refitted; no inference threshold was changed.",
             capture.output(table(primary$prefilter_status))), file.path(out,"audit_checks.txt"))
print(primary[!primary$retained_for_DE,])
cat("PASS prefilter reconstruction; no model refitting\n")
