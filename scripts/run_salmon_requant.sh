#!/usr/bin/env bash
set -euo pipefail

ORIGINAL=${ORIGINAL:-/path/to/rsoamp}
REVISION=${REVISION:-${ORIGINAL}/revision_R1_20260902}
SALMON=${SALMON:-/path/to/user-home/tools/miniconda3/envs/salmon_env/bin/salmon}
GFFREAD=${GFFREAD:-/path/to/user-home/tools/miniconda3/envs/starsolo/bin/gffread}
THREADS_PER_JOB=${THREADS_PER_JOB:-24}
MAX_JOBS=${MAX_JOBS:-2}

RAW=${REVISION}/04_unmapped_reads/raw_fastq
OUT=${REVISION}/06_expression_reanalysis
TRANSCRIPTS=${OUT}/rso_full_transcripts_from_gff.fa
DECOYS=${OUT}/rso_genome_decoys.txt
GENTROME=${OUT}/rso_transcripts_plus_genome_gentrome.fa
INDEX=${OUT}/salmon_index_decoy_aware_transcripts
QUANT=${OUT}/salmon_quant
LOG=${OUT}/logs
SAMPLES=${OUT}/rnaseq_sample_sheet.tsv
mkdir -p "${OUT}" "${QUANT}" "${LOG}" "${REVISION}/08_reproducibility/commands"

printf '%s\t%s\t%s\n' \
  run condition replicate \
  SRR27540881 CK 1 \
  SRR27540880 CK 2 \
  SRR27540879 CK 3 \
  SRR27540878 S200 1 \
  SRR27540877 S200 2 \
  SRR27540876 S200 3 \
  SRR27540875 S400 1 \
  SRR27540884 S400 2 \
  SRR27540883 S400 3 \
  > "${SAMPLES}"

pass_count=$(awk -F '\t' 'NR > 1 && $6 == "PASS" {n++} END {print n+0}' "${REVISION}/00_audit/rna_fastq_download_status.tsv")
[[ "${pass_count}" -eq 9 ]] || { echo "RNA FASTQ integrity gate is not PASS" >&2; exit 1; }

if [[ ! -s "${TRANSCRIPTS}" ]]; then
  rm -f "${TRANSCRIPTS}.tmp"
  "${GFFREAD}" "${ORIGINAL}/data/genome.gff" -g "${ORIGINAL}/data/genome.fa" \
    -w "${TRANSCRIPTS}.tmp" \
    > "${LOG}/gffread_transcripts.stdout.log" 2> "${LOG}/gffread_transcripts.stderr.log"
  [[ $(grep -c '^>' "${TRANSCRIPTS}.tmp") -eq 21791 ]] || {
    echo "Expected 21,791 reconstructed transcripts" >&2
    exit 1
  }
  mv "${TRANSCRIPTS}.tmp" "${TRANSCRIPTS}"
fi

if [[ ! -s "${DECOYS}" ]]; then
  grep '^>' "${ORIGINAL}/data/genome.fa" | cut -d ' ' -f 1 | sed 's/^>//' > "${DECOYS}"
fi
if [[ ! -s "${GENTROME}" ]]; then
  rm -f "${GENTROME}.tmp"
  cat "${TRANSCRIPTS}" "${ORIGINAL}/data/genome.fa" > "${GENTROME}.tmp"
  mv "${GENTROME}.tmp" "${GENTROME}"
fi

if [[ ! -s "${INDEX}/versionInfo.json" ]]; then
  rm -rf "${INDEX}"
  "${SALMON}" index -t "${GENTROME}" -d "${DECOYS}" -i "${INDEX}" -p 32 -k 31 --keepDuplicates \
    > "${LOG}/salmon_index.stdout.log" 2> "${LOG}/salmon_index.stderr.log"
fi
[[ -s "${INDEX}/versionInfo.json" && -s "${INDEX}/complete_ref_lens.bin" ]] || {
  echo "Salmon index is incomplete" >&2
  exit 1
}

quant_one() {
  local run=$1 condition=$2 replicate=$3
  local outdir=${QUANT}/${run}
  local r1=${RAW}/${run}_1.fastq.gz
  local r2=${RAW}/${run}_2.fastq.gz
  if [[ -s "${outdir}/QUANT_COMPLETE.PASS" ]]; then
    echo "SKIP ${run}; completion marker exists"
    return 0
  fi
  rm -rf "${outdir}"
  mkdir -p "${outdir}"
  cat > "${REVISION}/08_reproducibility/commands/${run}.salmon.command.txt" <<EOF
${SALMON} quant -i ${INDEX} -l A -1 ${r1} -2 ${r2} -p ${THREADS_PER_JOB} --validateMappings --seqBias --gcBias --writeUnmappedNames --dumpEq -o ${outdir}
EOF
  echo "[$(date --iso-8601=seconds)] START Salmon ${run} ${condition} replicate ${replicate}"
  "${SALMON}" quant -i "${INDEX}" -l A -1 "${r1}" -2 "${r2}" \
    -p "${THREADS_PER_JOB}" --validateMappings --seqBias --gcBias \
    --writeUnmappedNames --dumpEq -o "${outdir}" \
    > "${LOG}/${run}.salmon.stdout.log" 2> "${LOG}/${run}.salmon.stderr.log"
  [[ -s "${outdir}/quant.sf" && -s "${outdir}/aux_info/meta_info.json" \
      && -s "${outdir}/aux_info/unmapped_names.txt" ]] || {
    echo "Salmon output incomplete for ${run}" >&2
    return 1
  }
  date --iso-8601=seconds > "${outdir}/QUANT_COMPLETE.PASS"
  echo "[$(date --iso-8601=seconds)] PASS Salmon ${run}"
}

export -f quant_one
export REVISION SALMON THREADS_PER_JOB RAW INDEX QUANT LOG
while IFS=$'\t' read -r run condition replicate; do
  quant_one "${run}" "${condition}" "${replicate}" &
  while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do wait -n; done
done < <(tail -n +2 "${SAMPLES}")
wait

complete_count=$(find "${QUANT}" -name QUANT_COMPLETE.PASS | wc -l)
[[ "${complete_count}" -eq 9 ]] || { echo "Only ${complete_count}/9 Salmon runs completed" >&2; exit 1; }

python3 "${OUT}/summarize_salmon_quant.py"
sha256sum "${TRANSCRIPTS}" "${DECOYS}" "${SAMPLES}" "${QUANT}"/*/quant.sf \
  > "${OUT}/salmon_inputs_outputs.sha256"
date --iso-8601=seconds > "${OUT}/SALMON_REQUANT_COMPLETE.PASS"
echo "PASS all nine Salmon quantifications"
