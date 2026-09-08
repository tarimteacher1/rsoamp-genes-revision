#!/usr/bin/env bash
set -euo pipefail

REVISION=${REVISION:-/path/to/rsoamp/revision_R1_20260902}
WORK=${REVISION}/04_unmapped_reads
RAW=${WORK}/raw_fastq
QUANT=${REVISION}/06_expression_reanalysis/salmon_quant
ALIGN=${WORK}/star_alignments
REF=${WORK}/reference_databases
OUT=${WORK}/classification
BATCH=${WORK}/classification_batch
LOG=${WORK}/logs
SEQKIT=${SEQKIT:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/seqkit}
FASTP=${FASTP:-/path/to/user-home/tools/miniconda3/envs/salmon_env/bin/fastp}
BBDUK=${BBDUK:-/path/to/user-home/tools/miniconda3/envs/busco/opt/bbmap-39.52-0/bbduk.sh}
PIGZ=${PIGZ:-/usr/bin/pigz}
JAVA_BIN=${JAVA_BIN:-/path/to/user-home/tools/miniconda3/envs/busco/bin}
THREADS=${THREADS:-16}
export PATH="${JAVA_BIN}:${PATH}"
mkdir -p "${OUT}" "${LOG}"

[[ -s "${REF}/REFERENCE_DATABASES_COMPLETE.PASS" ]] || { echo "Reference database gate not passed" >&2; exit 1; }
[[ -s "${WORK}/RNA_GENOME_ALIGNMENT_COMPLETE.PASS" ]] || { echo "STAR genome-alignment gate not passed" >&2; exit 1; }
[[ -s "${REVISION}/06_expression_reanalysis/SALMON_REQUANT_COMPLETE.PASS" ]] || { echo "Salmon re-quantification gate not passed" >&2; exit 1; }

count_r1() {
  "${SEQKIT}" stats -T "$1" | awk 'NR==2 {print $4}'
}

write_ids() {
  "${SEQKIT}" seq -n -i "$1"
}

mapfile -t RUNS < <(tail -n +2 "${WORK}/rnaseq_sample_sheet.tsv" | cut -f1)
[[ "${#RUNS[@]}" -eq 9 ]] || { echo "Expected 9 RNA-seq runs, found ${#RUNS[@]}" >&2; exit 1; }

# Phase 1: prepare mutually exclusive decoy/quality categories per sample.
for run in "${RUNS[@]}"; do
  sample=${OUT}/${run}
  marker=${sample}/PREPROCESS_COMPLETE.PASS
  if [[ -s "${marker}" ]]; then
    echo "SKIP ${run}; preprocessing marker exists"
    continue
  fi
  rm -rf "${sample}"
  mkdir -p "${sample}"
  r1=${RAW}/${run}_1.fastq.gz
  r2=${RAW}/${run}_2.fastq.gz
  salmon_names=${QUANT}/${run}/aux_info/unmapped_names.txt
  [[ -s "${r1}" && -s "${r2}" && -s "${salmon_names}" ]] || {
    echo "Missing preprocessing input for ${run}" >&2
    exit 1
  }

  echo "[$(date --iso-8601=seconds)] START unmapped preprocessing ${run}"
  awk 'BEGIN {OFS="\t"} {count[$2]++} END {for (flag in count) print flag,count[flag]}' "${salmon_names}" \
    | sort > "${sample}/salmon_unmapped_flags.tsv"
  awk '$2 != "d" {print $1}' "${salmon_names}" | sort -u > "${sample}/nondecoy_unassigned_fragment_ids.txt"
  awk 'BEGIN {OFS="\t"} $2 == "d" {print $1,"genome_decoy_assigned"}' "${salmon_names}" \
    > "${sample}/fragment_categories.decoy.tsv"

  "${SEQKIT}" grep -f "${sample}/nondecoy_unassigned_fragment_ids.txt" -j "${THREADS}" "${r1}" \
    -o "${sample}/nondecoy_unassigned_1.fastq.gz"
  "${SEQKIT}" grep -f "${sample}/nondecoy_unassigned_fragment_ids.txt" -j "${THREADS}" "${r2}" \
    -o "${sample}/nondecoy_unassigned_2.fastq.gz"
  [[ "$(count_r1 "${sample}/nondecoy_unassigned_1.fastq.gz")" -eq "$(wc -l < "${sample}/nondecoy_unassigned_fragment_ids.txt")" ]] || {
    echo "Salmon name extraction mismatch for ${run}" >&2
    exit 1
  }

  "${FASTP}" \
    -i "${sample}/nondecoy_unassigned_1.fastq.gz" \
    -I "${sample}/nondecoy_unassigned_2.fastq.gz" \
    -o "${sample}/quality_pass_1.fastq.gz" \
    -O "${sample}/quality_pass_2.fastq.gz" \
    --failed_out "${sample}/quality_or_complexity_failed.fastq.gz" \
    --disable_adapter_trimming --disable_trim_poly_g \
    --qualified_quality_phred 15 --unqualified_percent_limit 40 \
    --length_required 50 --low_complexity_filter --complexity_threshold 30 \
    --thread "${THREADS}" \
    --json "${sample}/fastp_unmapped.json" --html "${sample}/fastp_unmapped.html" \
    > "${LOG}/${run}.unmapped_fastp.stdout.log" 2> "${LOG}/${run}.unmapped_fastp.stderr.log"
  write_ids "${sample}/quality_pass_1.fastq.gz" | sort -u > "${sample}/fastp_pass_fragment_ids.txt"
  write_ids "${sample}/quality_or_complexity_failed.fastq.gz" | sort -u > "${sample}/fastp_failed_fragment_ids.txt"
  cat "${sample}/fastp_pass_fragment_ids.txt" "${sample}/fastp_failed_fragment_ids.txt" \
    | sort -u > "${sample}/fastp_accounted_fragment_ids.txt"
  comm -23 "${sample}/nondecoy_unassigned_fragment_ids.txt" "${sample}/fastp_accounted_fragment_ids.txt" \
    > "${sample}/fastp_filtered_without_failed_out_ids.txt"
  cat "${sample}/fastp_failed_fragment_ids.txt" "${sample}/fastp_filtered_without_failed_out_ids.txt" \
    | sort -u | awk 'BEGIN {OFS="\t"} {print $1,"low_quality_or_low_complexity"}' \
    > "${sample}/fragment_categories.quality.tsv"
  date --iso-8601=seconds > "${marker}"
  echo "[$(date --iso-8601=seconds)] PASS unmapped preprocessing ${run}"
done

# Phase 2: concatenate valid gzip members, load each large reference once, then
# split mutually exclusive assignments back to samples by SRR-prefixed read ID.
if [[ ! -s "${BATCH}/BATCH_CLASSIFICATION_COMPLETE.PASS" ]]; then
  rm -rf "${BATCH}"
  mkdir -p "${BATCH}"
  # Recompress to one gzip stream per mate. BBMap's threaded gzip reader does
  # not reliably traverse directly concatenated gzip members.
  for mate in 1 2; do
    {
      for run in "${RUNS[@]}"; do
        "${PIGZ}" -dc "${OUT}/${run}/quality_pass_${mate}.fastq.gz"
      done
    } | "${PIGZ}" -p "${THREADS}" -c > "${BATCH}/quality_pass_${mate}.fastq.gz"
    "${PIGZ}" -t "${BATCH}/quality_pass_${mate}.fastq.gz"
  done

  "${BBDUK}" -Xmx64g \
    in1="${BATCH}/quality_pass_1.fastq.gz" in2="${BATCH}/quality_pass_2.fastq.gz" \
    ref="${REF}/SILVA_138.2_SSU_LSU_NR99.fasta" \
    outm1="${BATCH}/rrna_1.fastq.gz" outm2="${BATCH}/rrna_2.fastq.gz" \
    outu1="${BATCH}/nonrrna_1.fastq.gz" outu2="${BATCH}/nonrrna_2.fastq.gz" \
    k=27 hdist=0 removeifeitherbad=t overwrite=t threads="${THREADS}" \
    stats="${BATCH}/bbduk_rrna.stats.txt" \
    > "${LOG}/batch.bbduk_rrna.stdout.log" 2> "${LOG}/batch.bbduk_rrna.stderr.log"
  write_ids "${BATCH}/rrna_1.fastq.gz" \
    | awk 'BEGIN {OFS="\t"} {print $1,"rRNA_derived"}' \
    > "${BATCH}/fragment_categories.rrna.tsv"

  "${BBDUK}" -Xmx16g \
    in1="${BATCH}/nonrrna_1.fastq.gz" in2="${BATCH}/nonrrna_2.fastq.gz" \
    ref="${REF}/Reaumuria_songarica_chloroplast_NC_041273.fasta" \
    outm1="${BATCH}/plastid_1.fastq.gz" outm2="${BATCH}/plastid_2.fastq.gz" \
    outu1="${BATCH}/nonplastid_1.fastq.gz" outu2="${BATCH}/nonplastid_2.fastq.gz" \
    k=31 hdist=1 removeifeitherbad=t overwrite=t threads="${THREADS}" \
    stats="${BATCH}/bbduk_plastid.stats.txt" \
    > "${LOG}/batch.bbduk_plastid.stdout.log" 2> "${LOG}/batch.bbduk_plastid.stderr.log"
  write_ids "${BATCH}/plastid_1.fastq.gz" \
    | awk 'BEGIN {OFS="\t"} {print $1,"plastid_derived"}' \
    > "${BATCH}/fragment_categories.plastid.tsv"

  "${BBDUK}" -Xmx16g \
    in1="${BATCH}/nonplastid_1.fastq.gz" in2="${BATCH}/nonplastid_2.fastq.gz" \
    ref="${REF}/Myricaria_laxiflora_mitochondrion_MW971331_proxy.fasta" \
    outm1="${BATCH}/mitochondrial_proxy_1.fastq.gz" outm2="${BATCH}/mitochondrial_proxy_2.fastq.gz" \
    outu1="${BATCH}/residual_1.fastq.gz" outu2="${BATCH}/residual_2.fastq.gz" \
    k=31 hdist=1 removeifeitherbad=t overwrite=t threads="${THREADS}" \
    stats="${BATCH}/bbduk_mitochondrial_proxy.stats.txt" \
    > "${LOG}/batch.bbduk_mito.stdout.log" 2> "${LOG}/batch.bbduk_mito.stderr.log"
  write_ids "${BATCH}/mitochondrial_proxy_1.fastq.gz" \
    | awk 'BEGIN {OFS="\t"} {print $1,"mitochondrial_like_same_family_proxy"}' \
    > "${BATCH}/fragment_categories.mito.tsv"

  : > "${BATCH}/star_genome_unmapped_fragment_ids.txt"
  for run in "${RUNS[@]}"; do
    write_ids "${ALIGN}/${run}/${run}.Unmapped.out.mate1.gz" >> "${BATCH}/star_genome_unmapped_fragment_ids.txt"
  done
  sort -u -o "${BATCH}/star_genome_unmapped_fragment_ids.txt" "${BATCH}/star_genome_unmapped_fragment_ids.txt"
  "${SEQKIT}" grep -f "${BATCH}/star_genome_unmapped_fragment_ids.txt" -j "${THREADS}" \
    "${BATCH}/residual_1.fastq.gz" -o "${BATCH}/residual_genome_unmapped_1.fastq.gz"
  "${SEQKIT}" grep -f "${BATCH}/star_genome_unmapped_fragment_ids.txt" -j "${THREADS}" \
    "${BATCH}/residual_2.fastq.gz" -o "${BATCH}/residual_genome_unmapped_2.fastq.gz"
  "${SEQKIT}" grep -v -f "${BATCH}/star_genome_unmapped_fragment_ids.txt" -j "${THREADS}" \
    "${BATCH}/residual_1.fastq.gz" -o "${BATCH}/residual_genome_mapped_1.fastq.gz"
  write_ids "${BATCH}/residual_genome_unmapped_1.fastq.gz" \
    | awk 'BEGIN {OFS="\t"} {print $1,"genome_unmapped_residual"}' \
    > "${BATCH}/fragment_categories.genome_unmapped.tsv"
  write_ids "${BATCH}/residual_genome_mapped_1.fastq.gz" \
    | awk 'BEGIN {OFS="\t"} {print $1,"genome_mapped_but_transcriptome_unassigned"}' \
    > "${BATCH}/fragment_categories.genome_mapped.tsv"

  {
    for run in "${RUNS[@]}"; do
      cat "${OUT}/${run}/fragment_categories.decoy.tsv" "${OUT}/${run}/fragment_categories.quality.tsv"
    done
    cat "${BATCH}"/fragment_categories.rrna.tsv \
        "${BATCH}"/fragment_categories.plastid.tsv \
        "${BATCH}"/fragment_categories.mito.tsv \
        "${BATCH}"/fragment_categories.genome_mapped.tsv \
        "${BATCH}"/fragment_categories.genome_unmapped.tsv
  } | sort -k1,1 | gzip -c > "${BATCH}/fragment_category_assignments.tsv.gz"
  date --iso-8601=seconds > "${BATCH}/BATCH_CLASSIFICATION_COMPLETE.PASS"
fi

python3 "${WORK}/split_unmapped_batch_by_run.py" \
  --assignments "${BATCH}/fragment_category_assignments.tsv.gz" \
  --read1 "${BATCH}/residual_genome_unmapped_1.fastq.gz" \
  --read2 "${BATCH}/residual_genome_unmapped_2.fastq.gz" \
  --output-root "${OUT}" --runs "${RUNS[@]}" \
  --counts "${WORK}/per_sample_partition_integrity.tsv"

for run in "${RUNS[@]}"; do
  sample=${OUT}/${run}
  salmon_names=${QUANT}/${run}/aux_info/unmapped_names.txt
  assignments=$(awk -F '\t' -v run="${run}" '$1 == run {print $2}' "${WORK}/per_sample_partition_integrity.tsv")
  salmon_unassigned=$(wc -l < "${salmon_names}")
  [[ "${assignments}" -eq "${salmon_unassigned}" ]] || {
    echo "Category partition mismatch for ${run}: assignments=${assignments} salmon=${salmon_unassigned}" >&2
    exit 1
  }
  date --iso-8601=seconds > "${sample}/CLASSIFICATION_COMPLETE.PASS"
  echo "[$(date --iso-8601=seconds)] PASS unmapped classification ${run}"
done

python3 "${WORK}/summarize_unmapped_audit.py" \
  --classification-root "${OUT}" --quant-root "${QUANT}" --alignment-root "${ALIGN}" --output-dir "${WORK}"
date --iso-8601=seconds > "${WORK}/UNMAPPED_READ_CLASSIFICATION_COMPLETE.PASS"
echo "PASS all nine unmapped-read classifications"
