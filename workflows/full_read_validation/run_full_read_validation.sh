#!/usr/bin/env bash
set -euo pipefail

PROJECT=${PROJECT:-/path/to/rsoamp}
REVISION=${REVISION:-${PROJECT}/revision_R1_20260902}
SOURCE_RESULTS=${SOURCE_RESULTS:-${PROJECT}/r1_m6_mapped_rescue_20260905/run_retry1/results}
RUN_ROOT=${RUN_ROOT:-${PROJECT}/r1_m6_full_read_validation_20260905}
PREP=${PREP:-${RUN_ROOT}/prep}
OUT=${OUT:-${RUN_ROOT}/results}

STAR=${STAR:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/STAR}
SAMTOOLS=${SAMTOOLS:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/samtools}
STRINGTIE=${STRINGTIE:-/path/to/user-home/tools/miniconda3/envs/genome_annot/bin/stringtie}
GFFCOMPARE=${GFFCOMPARE:-/path/to/user-home/tools/miniconda3/envs/genome_annot/bin/gffcompare}
GFFREAD=${GFFREAD:-/path/to/user-home/tools/miniconda3/envs/commom_tools2/bin/gffread}
GETORF=${GETORF:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/getorf}
DIAMOND=${DIAMOND:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/diamond}
HMMSEARCH=${HMMSEARCH:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/hmmsearch}
TMBED_PY=${TMBED_PY:-${REVISION}/08_reproducibility/envs/tmbed/bin/python}

RAW=${REVISION}/04_unmapped_reads/raw_fastq
INDEX=${REVISION}/04_unmapped_reads/star_index_sjdb149
SAMPLES=${REVISION}/04_unmapped_reads/rnaseq_sample_sheet.tsv
REFERENCE_GTF=${REVISION}/04_unmapped_reads/reference_annotation.gtf
GENOME=${PROJECT}/data/genome.fa
WINDOWS=${PREP}/candidate_windows_10kb.bed
THREADS_STAR=${THREADS_STAR:-6}
THREADS_AUX=${THREADS_AUX:-1}

mkdir -p "${OUT}/align" "${OUT}/assembly" "${OUT}/screen" "${OUT}/logs" "${OUT}/commands"
exec > >(tee -a "${OUT}/logs/full_validation.master.log") 2>&1

required=(
  "${SOURCE_RESULTS}/candidate_evidence.tsv"
  "${SOURCE_RESULTS}/screen/references.dmnd"
  "${SOURCE_RESULTS}/screen/labeled_references.faa"
  "${PREP}/intergenic_candidate_models.gtf"
  "${PREP}/priority_start_stop_suborfs.faa"
  "${WINDOWS}"
  "${REFERENCE_GTF}"
  "${GENOME}"
  "${SAMPLES}"
)
for path in "${required[@]}"; do
  [[ -s "${path}" ]] || { echo "Missing required input: ${path}" >&2; exit 1; }
done
for path in Genome SA SAindex genomeParameters.txt; do
  [[ -s "${INDEX}/${path}" ]] || { echo "Incomplete STAR index: ${INDEX}/${path}" >&2; exit 1; }
done

pass_count=$(awk -F '\t' 'NR > 1 && $6 == "PASS" {n++} END {print n+0}' "${REVISION}/00_audit/rna_fastq_download_status.tsv")
[[ "${pass_count}" -eq 9 ]] || { echo "RNA FASTQ integrity gate is not 9/9 PASS" >&2; exit 1; }

sha256sum \
  "${SOURCE_RESULTS}/candidate_evidence.tsv" \
  "${SOURCE_RESULTS}/screen/comparison.annotated.gtf" \
  "${PREP}/intergenic_candidate_models.gtf" \
  "${PREP}/priority_start_stop_suborfs.faa" \
  "${WINDOWS}" \
  "${REFERENCE_GTF}" \
  > "${OUT}/input_identities.sha256"

align_one() {
  local run=$1 treatment=$2 replicate=$3
  local dir=${OUT}/align/${run}
  local prefix=${dir}/${run}.
  local bam=${dir}/${run}.candidate_windows.bam
  local tmp=${bam}.tmp
  local r1=${RAW}/${run}_1.fastq.gz
  local r2=${RAW}/${run}_2.fastq.gz
  mkdir -p "${dir}"
  if [[ -s "${dir}/FULL_READ_LOCAL_ALIGNMENT.PASS" ]]; then
    echo "[$(date --iso-8601=seconds)] SKIP ${run}; PASS marker exists"
    return 0
  fi
  [[ -s "${r1}" && -s "${r2}" ]] || { echo "Missing FASTQ pair for ${run}" >&2; return 1; }
  rm -f "${tmp}" "${tmp}.bai" "${bam}" "${bam}.bai"
  cat > "${OUT}/commands/${run}.STAR_stream.command.txt" <<EOF
${STAR} --runThreadN ${THREADS_STAR} --genomeDir ${INDEX} --genomeLoad NoSharedMemory --readFilesIn ${r1} ${r2} --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix ${prefix} --outSAMtype BAM Unsorted --outStd BAM_Unsorted --outSAMunmapped Within KeepPairs --outSAMattributes NH HI AS nM XS --outSAMstrandField intronMotif --outSAMattrRGline ID:${run} SM:${run} | ${SAMTOOLS} view -@ ${THREADS_AUX} -bh -L ${WINDOWS} - | ${SAMTOOLS} sort -@ ${THREADS_AUX} -o ${bam} -
EOF
  echo "[$(date --iso-8601=seconds)] START full-read local alignment ${run} ${treatment} replicate ${replicate}"
  set +e
  "${STAR}" \
    --runThreadN "${THREADS_STAR}" \
    --genomeDir "${INDEX}" \
    --genomeLoad NoSharedMemory \
    --readFilesIn "${r1}" "${r2}" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFileNamePrefix "${prefix}" \
    --outSAMtype BAM Unsorted \
    --outStd BAM_Unsorted \
    --outSAMunmapped Within KeepPairs \
    --outSAMattributes NH HI AS nM XS \
    --outSAMstrandField intronMotif \
    --outSAMattrRGline "ID:${run}" "SM:${run}" \
    2> "${dir}/${run}.STAR.stderr.log" \
    | "${SAMTOOLS}" view -@ "${THREADS_AUX}" -bh -L "${WINDOWS}" - \
    | "${SAMTOOLS}" sort -@ "${THREADS_AUX}" -o "${tmp}" -
  statuses=("${PIPESTATUS[@]}")
  set -e
  if [[ "${statuses[0]}" -ne 0 || "${statuses[1]}" -ne 0 || "${statuses[2]}" -ne 0 ]]; then
    echo "Pipeline failure for ${run}: STAR=${statuses[0]} view=${statuses[1]} sort=${statuses[2]}" >&2
    return 1
  fi
  mv "${tmp}" "${bam}"
  "${SAMTOOLS}" quickcheck -v "${bam}"
  "${SAMTOOLS}" index -@ "${THREADS_AUX}" "${bam}"
  "${SAMTOOLS}" flagstat -O json "${bam}" > "${dir}/${run}.candidate_windows.flagstat.json"
  "${SAMTOOLS}" bedcov "${WINDOWS}" "${bam}" > "${dir}/${run}.candidate_windows.bedcov.tsv"

  observed=$(awk -F '|' '/Number of input reads/ {gsub(/[[:space:]]/, "", $2); print $2}' "${prefix}Log.final.out")
  expected=$(awk -F '\t' -v run="${run}" '$1 == run {print $2}' "${REVISION}/00_audit/rna_fastq_download_status.tsv")
  [[ -n "${observed}" && "${observed}" == "${expected}" ]] || {
    echo "Raw-read count mismatch for ${run}: STAR=${observed:-missing}, audit=${expected:-missing}" >&2
    return 1
  }
  printf 'run\ttreatment\treplicate\texpected_input_fragments\tstar_input_fragments\tretained_alignment_records\n%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${run}" "${treatment}" "${replicate}" "${expected}" "${observed}" \
    "$("${SAMTOOLS}" view -c "${bam}")" > "${dir}/${run}.integrity.tsv"

  "${STRINGTIE}" "${bam}" -G "${REFERENCE_GTF}" --rf -p 2 -m 90 -c 1 -s 2 -f 0.01 -j 2 -a 10 -M 0.5 \
    -l "M6FULL_${run}" -o "${OUT}/assembly/${run}.gtf" \
    > "${OUT}/logs/${run}.StringTie.stdout.log" 2> "${OUT}/logs/${run}.StringTie.stderr.log"
  [[ -s "${OUT}/assembly/${run}.gtf" ]] || { echo "Empty StringTie GTF for ${run}" >&2; return 1; }
  date --iso-8601=seconds > "${dir}/FULL_READ_LOCAL_ALIGNMENT.PASS"
  echo "[$(date --iso-8601=seconds)] PASS full-read local alignment ${run}"
}

while IFS=$'\t' read -r run treatment replicate; do
  align_one "${run}" "${treatment}" "${replicate}"
done < <(tail -n +2 "${SAMPLES}")

find "${OUT}/assembly" -maxdepth 1 -type f -name '*.gtf' -print | sort > "${OUT}/screen/assemblies.txt"
[[ "$(wc -l < "${OUT}/screen/assemblies.txt")" -eq 9 ]] || { echo "Expected nine local assemblies" >&2; exit 1; }
"${STRINGTIE}" --merge -p 4 -m 90 -f 0.01 -l M6FULLMERGE -o "${OUT}/screen/full_read_merged.gtf" "${OUT}/screen/assemblies.txt"
"${GFFCOMPARE}" -r "${PREP}/intergenic_candidate_models.gtf" -o "${OUT}/screen/vs_subset_candidate" "${OUT}/screen/full_read_merged.gtf"
"${GFFCOMPARE}" -r "${REFERENCE_GTF}" -o "${OUT}/screen/vs_reference" "${OUT}/screen/full_read_merged.gtf"
"${GFFREAD}" "${OUT}/screen/full_read_merged.gtf" -g "${GENOME}" -w "${OUT}/screen/full_read_transcripts.fna"
PATH="$(dirname "${GETORF}"):${PATH}" "${GETORF}" -sequence "${OUT}/screen/full_read_transcripts.fna" \
  -outseq "${OUT}/screen/full_read_stop_delimited_orfs.faa" -find 0 -minsize 90 -reverse N -table 0 -auto
"${DIAMOND}" blastp --query "${OUT}/screen/full_read_stop_delimited_orfs.faa" \
  --db "${SOURCE_RESULTS}/screen/references.dmnd" --out "${OUT}/screen/full_read_orf_reference_hits.tsv" \
  --sensitive --evalue 0.001 --id 25 --max-target-seqs 0 --threads 6 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore
"${HMMSEARCH}" --cpu 6 --cut_ga --noali --domtblout "${OUT}/screen/full_read_orf_family_domains.domtbl" \
  "${PROJECT}/data/hmm/AMP_models.hmm" "${OUT}/screen/full_read_stop_delimited_orfs.faa" \
  > "${OUT}/logs/full_read_hmmsearch.log"

"${DIAMOND}" blastp --query "${PREP}/priority_start_stop_suborfs.faa" \
  --db "${SOURCE_RESULTS}/screen/references.dmnd" --out "${OUT}/screen/priority_vs_labeled_AMP_refs.tsv" \
  --sensitive --evalue 0.001 --id 20 --max-target-seqs 0 --threads 6 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore
"${DIAMOND}" blastp --query "${PREP}/priority_start_stop_suborfs.faa" \
  --db "${PROJECT}/intermediate/mcscanx/rso.dmnd" --out "${OUT}/screen/priority_vs_rso_proteome.tsv" \
  --sensitive --evalue 0.001 --id 20 --max-target-seqs 25 --threads 6 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore
"${HMMSEARCH}" --cpu 6 --cut_ga --noali --domtblout "${OUT}/screen/priority_family_domains.domtbl" \
  "${PROJECT}/data/hmm/AMP_models.hmm" "${PREP}/priority_start_stop_suborfs.faa" \
  > "${OUT}/logs/priority_hmmsearch.log"

if [[ -x "${TMBED_PY}" && -d "${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc" ]]; then
  export HF_HOME="${REVISION}/08_reproducibility/models/huggingface_cache"
  export TRANSFORMERS_OFFLINE=1
  export HF_HUB_OFFLINE=1
  "${TMBED_PY}" -m tmbed predict \
    --fasta "${PREP}/priority_start_stop_suborfs.faa" \
    --predictions "${OUT}/screen/priority_start_stop_suborfs.tmbed.pred" \
    --out-format 0 --no-use-gpu --threads 4 \
    --model-dir "${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc" \
    > "${OUT}/logs/priority_tmbed.log" 2>&1
else
  echo "TMbed runtime unavailable" > "${OUT}/screen/TMBED_NOT_RUN.txt"
fi

date --iso-8601=seconds > "${OUT}/FULL_READ_LOCAL_VALIDATION_RAW_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS full-read local validation raw workflow"
