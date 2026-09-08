#!/usr/bin/env bash
set -euo pipefail

ORIGINAL=${ORIGINAL:-/path/to/rsoamp}
REVISION=${REVISION:-${ORIGINAL}/revision_R1_20260902}
STAR=${STAR:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/STAR}
SAMTOOLS=${SAMTOOLS:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/samtools}
THREADS_PER_JOB=${THREADS_PER_JOB:-32}
MAX_JOBS=${MAX_JOBS:-2}

WORK=${REVISION}/04_unmapped_reads
RAW=${WORK}/raw_fastq
INDEX=${WORK}/star_index_sjdb149
ALIGN=${WORK}/star_alignments
TARGET=${WORK}/backfill_targets_10kb.bed
RUNS=${WORK}/rnaseq_sample_sheet.tsv
LOGDIR=${WORK}/logs
mkdir -p "${ALIGN}" "${LOGDIR}" "${REVISION}/08_reproducibility/commands"

[[ -s "${INDEX}/INDEX_COMPLETE.PASS" ]] || {
  echo "STAR index is incomplete: ${INDEX}" >&2
  exit 1
}

python3 - "${REVISION}/03_annotation_rescue/backfill_locus_audit_preliminary.tsv" "${TARGET}" <<'PY'
import csv
import sys
from pathlib import Path

source, target = map(Path, sys.argv[1:])
with source.open(encoding="utf-8") as handle, target.open("w", encoding="utf-8") as out:
    for row in csv.DictReader(handle, delimiter="\t"):
        start = max(0, int(row["coding_start"]) - 10001)
        end = int(row["coding_end"]) + 10000
        out.write(f'{row["historical_locus"].split(":", 1)[0]}\t{start}\t{end}\t{row["member_id"]}\n')
PY

printf '%s\t%s\t%s\n' \
  run treatment replicate \
  SRR27540881 CK 1 \
  SRR27540880 CK 2 \
  SRR27540879 CK 3 \
  SRR27540878 S200 1 \
  SRR27540877 S200 2 \
  SRR27540876 S200 3 \
  SRR27540875 S400 1 \
  SRR27540884 S400 2 \
  SRR27540883 S400 3 \
  > "${RUNS}"

pass_count=$(awk -F '\t' 'NR > 1 && $6 == "PASS" {n++} END {print n+0}' "${REVISION}/00_audit/rna_fastq_download_status.tsv")
[[ "${pass_count}" -eq 9 ]] || {
  echo "Expected nine PASS rows in RNA FASTQ integrity table; observed ${pass_count}" >&2
  exit 1
}

align_one() {
  local run=$1 treatment=$2 replicate=$3
  local outdir=${ALIGN}/${run}
  local prefix=${outdir}/${run}.
  local log=${LOGDIR}/${run}.STAR.log
  local r1=${RAW}/${run}_1.fastq.gz
  local r2=${RAW}/${run}_2.fastq.gz
  mkdir -p "${outdir}"

  if [[ -s "${outdir}/ALIGNMENT_COMPLETE.PASS" ]]; then
    echo "[$(date --iso-8601=seconds)] SKIP ${run}; completion marker exists"
    return 0
  fi
  [[ -s "${r1}" && -s "${r2}" ]] || {
    echo "Missing FASTQ pair for ${run}" >&2
    return 1
  }

  rm -rf "${outdir}"
  mkdir -p "${outdir}"
  echo "[$(date --iso-8601=seconds)] START ${run} ${treatment} replicate ${replicate}"
  cat > "${REVISION}/08_reproducibility/commands/${run}.STAR.command.txt" <<EOF
${STAR} --runThreadN ${THREADS_PER_JOB} --genomeDir ${INDEX} --genomeLoad NoSharedMemory --readFilesIn ${r1} ${r2} --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ${prefix} --outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs --outReadsUnmapped Fastx --outSAMattributes NH HI AS nM XS --outSAMstrandField intronMotif --quantMode GeneCounts
EOF

  "${STAR}" \
    --runThreadN "${THREADS_PER_JOB}" \
    --genomeDir "${INDEX}" \
    --genomeLoad NoSharedMemory \
    --readFilesIn "${r1}" "${r2}" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFileNamePrefix "${prefix}" \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within KeepPairs \
    --outReadsUnmapped Fastx \
    --outSAMattributes NH HI AS nM XS \
    --outSAMstrandField intronMotif \
    --quantMode GeneCounts 2>&1 | tee "${log}"

  [[ -s "${prefix}Aligned.out.bam" && -s "${prefix}Log.final.out" && -s "${prefix}SJ.out.tab" ]] || {
    echo "Incomplete STAR outputs for ${run}" >&2
    return 1
  }

  "${SAMTOOLS}" view -@ 4 -h -L "${TARGET}" "${prefix}Aligned.out.bam" \
    | "${SAMTOOLS}" sort -@ 4 -o "${outdir}/${run}.backfill_targets.bam" -
  "${SAMTOOLS}" index "${outdir}/${run}.backfill_targets.bam"
  "${SAMTOOLS}" bedcov "${TARGET}" "${outdir}/${run}.backfill_targets.bam" \
    > "${outdir}/${run}.backfill_targets.bedcov.tsv"
  "${SAMTOOLS}" coverage "${outdir}/${run}.backfill_targets.bam" \
    > "${outdir}/${run}.backfill_targets.contig_coverage.tsv"

  gzip -f "${prefix}Unmapped.out.mate1" "${prefix}Unmapped.out.mate2"
  rm -f "${prefix}Aligned.out.bam"

  {
    printf 'run\ttreatment\treplicate\ttarget_bam_alignments\tunmapped_mate1_reads\tunmapped_mate2_reads\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${run}" "${treatment}" "${replicate}" \
      "$("${SAMTOOLS}" view -c "${outdir}/${run}.backfill_targets.bam")" \
      "$(zcat "${prefix}Unmapped.out.mate1.gz" | awk 'END {print NR/4}')" \
      "$(zcat "${prefix}Unmapped.out.mate2.gz" | awk 'END {print NR/4}')"
  } > "${outdir}/${run}.alignment_integrity.tsv"

  date --iso-8601=seconds > "${outdir}/ALIGNMENT_COMPLETE.PASS"
  echo "[$(date --iso-8601=seconds)] PASS ${run}"
}

export -f align_one
export ORIGINAL REVISION STAR SAMTOOLS THREADS_PER_JOB RAW INDEX ALIGN TARGET LOGDIR

while IFS=$'\t' read -r run treatment replicate; do
  align_one "${run}" "${treatment}" "${replicate}" &
  while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do
    wait -n
  done
done < <(tail -n +2 "${RUNS}")
wait

{
  head -1 "${ALIGN}/SRR27540881/SRR27540881.alignment_integrity.tsv"
  find "${ALIGN}" -name '*.alignment_integrity.tsv' -print0 \
    | sort -z \
    | xargs -0 -n1 tail -n +2
} > "${WORK}/rna_genome_alignment_integrity.tsv"

complete_count=$(find "${ALIGN}" -name ALIGNMENT_COMPLETE.PASS | wc -l)
[[ "${complete_count}" -eq 9 ]] || {
  echo "Only ${complete_count}/9 STAR alignments have PASS markers" >&2
  exit 1
}
date --iso-8601=seconds > "${WORK}/RNA_GENOME_ALIGNMENT_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS all nine RNA-seq genome alignments"
