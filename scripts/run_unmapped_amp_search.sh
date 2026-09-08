#!/usr/bin/env bash
set -euo pipefail

REVISION=${REVISION:-/path/to/rsoamp/revision_R1_20260902}
WORK=${REVISION}/04_unmapped_reads
CLASS=${WORK}/classification
REF=${WORK}/reference_databases
DIAMOND=${DIAMOND:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/diamond}
THREADS=${THREADS:-24}

[[ -s "${WORK}/UNMAPPED_READ_CLASSIFICATION_COMPLETE.PASS" ]] || { echo "Unmapped classification gate not passed" >&2; exit 1; }
[[ -s "${REVISION}/03_annotation_rescue/FINAL_AMP_CATALOGUE_COMPLETE.PASS" ]] || { echo "Final AMP catalogue gate not passed" >&2; exit 1; }

python3 "${WORK}/build_unmapped_amp_reference.py" \
  --revision "${REVISION}" \
  --output "${REF}/targeted_AMP_reference.faa" \
  --manifest "${REF}/targeted_AMP_reference_manifest.tsv"
"${DIAMOND}" makedb --in "${REF}/targeted_AMP_reference.faa" --db "${REF}/targeted_AMP_reference" \
  > "${WORK}/logs/diamond_AMP_makedb.log" 2>&1

while IFS=$'\t' read -r run treatment replicate; do
  sample=${CLASS}/${run}
  query=${sample}/residual_genome_unmapped_reads.fasta.gz
  python3 "${WORK}/fastq_pairs_to_labeled_fasta.py" \
    --read1 "${sample}/residual_genome_unmapped_1.fastq.gz" \
    --read2 "${sample}/residual_genome_unmapped_2.fastq.gz" \
    --output "${query}"
  "${DIAMOND}" blastx \
    --query "${query}" --db "${REF}/targeted_AMP_reference" \
    --out "${sample}/residual_genome_unmapped_vs_AMP.diamond.tsv" \
    --more-sensitive --evalue 1e-5 --id 50 --query-cover 50 --subject-cover 30 \
    --max-target-seqs 25 --threads "${THREADS}" \
    --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore \
    > "${WORK}/logs/${run}.diamond_AMP.stdout.log" 2> "${WORK}/logs/${run}.diamond_AMP.stderr.log"
done < <(tail -n +2 "${WORK}/rnaseq_sample_sheet.tsv")

python3 "${WORK}/summarize_unmapped_amp_hits.py" \
  --classification-root "${CLASS}" \
  --output "${WORK}/unmapped_amp_hits.tsv" \
  --summary "${WORK}/unmapped_amp_search_summary.tsv"
date --iso-8601=seconds > "${WORK}/UNMAPPED_AMP_SEARCH_COMPLETE.PASS"
echo "PASS targeted AMP search"
