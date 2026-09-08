#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(apeglm)
})

original <- "/path/to/rsoamp"
revision <- file.path(original, "revision_R1_20260902")
out <- file.path(revision, "06_expression_reanalysis")
dir.create(file.path(out, "results"), recursive=TRUE, showWarnings=FALSE)

samples <- read.delim(file.path(out, "rnaseq_sample_sheet.tsv"), stringsAsFactors=FALSE)
samples$condition <- factor(samples$condition, levels=c("CK", "S200", "S400"))
rownames(samples) <- samples$run
files <- file.path(out, "salmon_quant", samples$run, "quant.sf")
names(files) <- samples$run
stopifnot(all(file.exists(files)))

gff <- readLines(file.path(original, "data/genome.gff"), warn=FALSE)
mrna <- strsplit(gff[grepl("\\tmRNA\\t", gff)], "\\t")
extract_attr <- function(text, key) {
  hit <- regmatches(text, regexpr(paste0("(?:^|;)", key, "=[^;]+"), text, perl=TRUE))
  sub(paste0("^(?:.*;)?", key, "="), "", hit)
}
tx2gene <- do.call(rbind, lapply(mrna, function(fields) {
  data.frame(TXNAME=extract_attr(fields[9], "ID"), GENEID=extract_attr(fields[9], "Parent"), stringsAsFactors=FALSE)
}))
tx2gene <- unique(tx2gene[nzchar(tx2gene$TXNAME) & nzchar(tx2gene$GENEID),])
stopifnot(!anyDuplicated(tx2gene$TXNAME))
write.table(tx2gene, file.path(out, "transcript_to_gene.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="no")
coldata <- samples[,c("condition", "replicate"),drop=FALSE]
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~condition)
pre_counts <- counts(dds)
keep <- rowSums(pre_counts) >= 10
prefilter <- data.frame(
  gene_id=rownames(pre_counts),
  tximport_estimated_count_sum=rowSums(txi$counts),
  deseq2_rounded_count_sum=rowSums(pre_counts),
  nonzero_count_samples=rowSums(pre_counts > 0),
  nonzero_tpm_samples=rowSums(txi$abundance > 0),
  max_TPM=apply(txi$abundance, 1, max),
  prefilter_threshold_total_count=10L,
  retained_for_DE=keep,
  prefilter_status=ifelse(keep, "RETAINED_TOTAL_COUNT_GE_10",
    ifelse(rowSums(pre_counts) == 0, "FILTERED_ALL_ZERO_COUNTS", "FILTERED_TOTAL_COUNT_LT_10"))
)
write.table(prefilter, file.path(out,"results/all_gene_prefilter_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene_id=rownames(pre_counts),pre_counts,check.names=FALSE),
  file.path(out,"results/all_gene_prefilter_counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
dds <- dds[keep,]
dds <- DESeq(dds)
saveRDS(dds, file.path(out, "results/deseq2_dataset.rds"))

raw_result <- function(a, b, tag) {
  result <- as.data.frame(results(dds, contrast=c("condition", a, b), alpha=0.05))
  result$gene_id <- rownames(result)
  result$contrast <- tag
  result
}
all_results <- rbind(
  raw_result("S200", "CK", "S200_vs_CK"),
  raw_result("S400", "CK", "S400_vs_CK"),
  raw_result("S400", "S200", "S400_vs_S200")
)

shrink_map <- list(S200_vs_CK="condition_S200_vs_CK", S400_vs_CK="condition_S400_vs_CK")
all_results$log2FC_apeglm <- NA_real_
for (tag in names(shrink_map)) {
  shrunk <- as.data.frame(lfcShrink(dds, coef=shrink_map[[tag]], type="apeglm", quiet=TRUE))
  indices <- which(all_results$contrast == tag)
  all_results$log2FC_apeglm[indices] <- shrunk$log2FoldChange[match(all_results$gene_id[indices], rownames(shrunk))]
}
all_results$DE_call <- !is.na(all_results$padj) & all_results$padj < 0.05 & abs(all_results$log2FoldChange) > 1
all_results <- all_results[,c("gene_id","contrast","baseMean","log2FoldChange","log2FC_apeglm","lfcSE","stat","pvalue","padj","DE_call")]
write.csv(all_results, file.path(out, "results/all_gene_DE_results.csv"), row.names=FALSE)

members <- read.csv(file.path(original, "tables/T1_AMP_members.csv"), stringsAsFactors=FALSE)
members <- members[members$source == "annotated",]
members$gene_id <- tx2gene$GENEID[match(members$protein_id, tx2gene$TXNAME)]
stopifnot(!any(is.na(members$gene_id)))
amp <- merge(members[,c("member_id","family","subclass","protein_id","gene_id")], all_results, by="gene_id", all.x=TRUE)
amp <- amp[order(amp$family, amp$member_id, amp$contrast),]
write.csv(amp, file.path(out, "results/annotated_AMP_DE_results_pre_catalogue_freeze.csv"), row.names=FALSE)

tpm <- txi$abundance
write.csv(data.frame(gene_id=rownames(tpm), tpm, check.names=FALSE), file.path(out, "results/gene_TPM.csv"), row.names=FALSE)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(data.frame(gene_id=rownames(normalized_counts), normalized_counts, check.names=FALSE), file.path(out, "results/gene_normalized_counts.csv"), row.names=FALSE)
vst_mat <- assay(vst(dds, blind=FALSE))
write.csv(data.frame(gene_id=rownames(vst_mat), vst_mat, check.names=FALSE), file.path(out, "results/gene_VST.csv"), row.names=FALSE)

pca <- plotPCA(vst(dds, blind=FALSE), intgroup="condition", returnData=TRUE)
pca$run <- rownames(pca)
write.csv(pca, file.path(out, "results/PCA_coordinates.csv"), row.names=FALSE)
write.csv(as.matrix(dist(t(vst_mat))), file.path(out, "results/sample_distance_matrix.csv"))
write.csv(cor(vst_mat, method="pearson"), file.path(out, "results/sample_pearson_correlation.csv"))

summary_rows <- do.call(rbind, lapply(unique(amp$contrast[!is.na(amp$contrast)]), function(tag) {
  subset <- amp[amp$contrast == tag & !is.na(amp$padj),]
  data.frame(
    contrast=tag,
    n_tested=nrow(subset),
    n_DE=sum(subset$DE_call),
    n_up=sum(subset$DE_call & subset$log2FoldChange > 1),
    n_down=sum(subset$DE_call & subset$log2FoldChange < -1)
  )
}))
write.table(summary_rows, file.path(out, "results/annotated_AMP_DE_summary_pre_catalogue_freeze.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

session <- capture.output(sessionInfo())
writeLines(session, file.path(out, "expression_R_sessionInfo.txt"))
cat("PASS tximport + DESeq2 expression reanalysis\n")
print(summary_rows)
