#!/usr/bin/env bash
set -euo pipefail

PROJECT=${PROJECT:-/path/to/rsoamp}
REVISION=${REVISION:-${PROJECT}/revision_R1_20260902}
SOURCE_RESULTS=${SOURCE_RESULTS:-${PROJECT}/r1_m6_mapped_rescue_20260905/run_retry1/results}
RUN_ROOT=${RUN_ROOT:-${PROJECT}/r1_m6_recursive_family_audit_20260905}
PREP=${PREP:-${RUN_ROOT}/prep}
OUT=${OUT:-${RUN_ROOT}/results}

DIAMOND=${DIAMOND:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/diamond}
HMMSEARCH=${HMMSEARCH:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/hmmsearch}
BLASTP=${BLASTP:-/path/to/user-home/tools/miniconda3/bin/blastp}
TMBED_PY=${TMBED_PY:-${REVISION}/08_reproducibility/envs/tmbed/bin/python}
PREDGPI=${PREDGPI:-${REVISION}/08_reproducibility/software/predgpi}
DEEPSIG=${DEEPSIG:-/path/to/user-home/tools/miniconda3/envs/deepsig_env/bin/deepsig}
DEEPSIG_ROOT=${DEEPSIG_ROOT:-/path/to/user-home/tools/miniconda3/envs/deepsig_env/share/deepsig}
NSLTP_LABELED_REFERENCE=${REVISION}/01_nslTP_curation/reference/combined_labeled_reference.faa
NSLTP_REFERENCE_SUMMARIZER=${REVISION}/05_comparative_genomics/summarize_labeled_reference_blast.py
QUERY=${PREP}/recursive_pf00234_start_stop.faa

mkdir -p "${OUT}/logs"
exec > >(tee -a "${OUT}/logs/recursive_family_audit.master.log") 2>&1
for path in "${QUERY}" "${PREP}/recursive_pf00234_manifest.tsv" "${NSLTP_LABELED_REFERENCE}"; do
  [[ -s "${path}" ]] || { echo "Missing required input: ${path}" >&2; exit 1; }
done

"${DIAMOND}" blastp --query "${QUERY}" \
  --db "${SOURCE_RESULTS}/screen/references.dmnd" \
  --out "${OUT}/recursive_vs_labeled_AMP_refs.tsv" \
  --sensitive --evalue 0.001 --id 20 --max-target-seqs 0 --threads 1 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore

"${DIAMOND}" blastp --query "${QUERY}" \
  --db "${PROJECT}/intermediate/mcscanx/rso.dmnd" \
  --out "${OUT}/recursive_vs_rso_proteome.tsv" \
  --sensitive --evalue 0.001 --id 20 --max-target-seqs 25 --threads 1 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore

"${HMMSEARCH}" --cpu 1 --cut_ga --noali \
  --domtblout "${OUT}/recursive_family_domains.domtbl" \
  "${PROJECT}/data/hmm/AMP_models.hmm" "${QUERY}" \
  > "${OUT}/logs/recursive_hmmsearch.log"

"${BLASTP}" -task blastp-short \
  -query "${QUERY}" -subject "${NSLTP_LABELED_REFERENCE}" \
  -evalue 10 -seg no -comp_based_stats 0 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
  -out "${OUT}/recursive_vs_nsltp_labeled_reference.tsv"
python3 "${NSLTP_REFERENCE_SUMMARIZER}" \
  --blast "${OUT}/recursive_vs_nsltp_labeled_reference.tsv" \
  --queries "${QUERY}" \
  --output "${OUT}/recursive_nsltp_reference_similarity.tsv"

if [[ -f "${PREDGPI}/predgpi.py" ]]; then
  export PREDGPI_HOME="${PREDGPI}"
  python3 "${PREDGPI}/predgpi.py" -f "${QUERY}" \
    -o "${OUT}/recursive_pf00234.predgpi.json" -m json \
    > "${OUT}/logs/recursive_predgpi.stdout.log" \
    2> "${OUT}/logs/recursive_predgpi.stderr.log"
else
  echo "PredGPI runtime unavailable" > "${OUT}/PREDGPI_NOT_RUN.txt"
fi

if [[ -x "${DEEPSIG}" && -d "${DEEPSIG_ROOT}" ]]; then
  export DEEPSIG_ROOT
  "${DEEPSIG}" -f "${QUERY}" -o "${OUT}/recursive_pf00234.deepsig.gff3" \
    -k euk -m gff3 -t 1 \
    > "${OUT}/logs/recursive_deepsig.stdout.log" \
    2> "${OUT}/logs/recursive_deepsig.stderr.log"
else
  echo "DeepSig runtime unavailable" > "${OUT}/DEEPSIG_NOT_RUN.txt"
fi

if [[ -x "${TMBED_PY}" && -d "${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc" ]]; then
  export HF_HOME="${REVISION}/08_reproducibility/models/huggingface_cache"
  export TRANSFORMERS_OFFLINE=1
  export HF_HUB_OFFLINE=1
  "${TMBED_PY}" -m tmbed predict \
    --fasta "${QUERY}" --predictions "${OUT}/recursive_pf00234.tmbed.pred" \
    --out-format 0 --no-use-gpu --threads 1 \
    --model-dir "${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc" \
    > "${OUT}/logs/recursive_tmbed.log" 2>&1
else
  echo "TMbed runtime unavailable" > "${OUT}/TMBED_NOT_RUN.txt"
fi

sha256sum "${QUERY}" "${PREP}/recursive_pf00234_manifest.tsv" \
  "${NSLTP_LABELED_REFERENCE}" > "${OUT}/input_identities.sha256"
date --iso-8601=seconds > "${OUT}/RECURSIVE_FAMILY_AUDIT_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS recursive family audit"
