#!/usr/bin/env bash
set -euo pipefail

PROJECT=${PROJECT:-/path/to/rsoamp}
REVISION=${REVISION:-${PROJECT}/revision_R1_20260902}
SOURCE_RESULTS=${SOURCE_RESULTS:-${PROJECT}/r1_m6_mapped_rescue_20260905/run_retry1/results}
RUN_ROOT=${RUN_ROOT:-${PROJECT}/r1_m6_overlap_model_audit_20260905}
PREP=${PREP:-${RUN_ROOT}/prep}
OUT=${OUT:-${RUN_ROOT}/results}
WORKFLOW=${WORKFLOW:-${RUN_ROOT}/workflow}

PYTHON=${PYTHON:-/path/to/user-home/tools/miniconda3/envs/genome_annot/bin/python}
DIAMOND=${DIAMOND:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/diamond}
HMMSEARCH=${HMMSEARCH:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/hmmsearch}
BLASTP=${BLASTP:-/path/to/user-home/tools/miniconda3/bin/blastp}
TMBED_PY=${TMBED_PY:-${REVISION}/08_reproducibility/envs/tmbed/bin/python}
PREDGPI=${PREDGPI:-${REVISION}/08_reproducibility/software/predgpi}
NSLTP_LABELED_REFERENCE=${REVISION}/01_nslTP_curation/reference/combined_labeled_reference.faa
NSLTP_REFERENCE_SUMMARIZER=${REVISION}/05_comparative_genomics/summarize_labeled_reference_blast.py
CATALOGUE_ALL=${REVISION}/03_annotation_rescue/final_amp_catalogue_all_audited_members.tsv
NSLTP_EVIDENCE=${REVISION}/01_nslTP_curation/nsLTP_evidence_matrix.tsv

mkdir -p "${PREP}" "${OUT}" "${OUT}/logs"
exec > >(tee -a "${OUT}/logs/overlap_family_audit.master.log") 2>&1

"${PYTHON}" "${WORKFLOW}/prepare_overlap_model_audit.py" \
  --source-results "${SOURCE_RESULTS}" \
  --catalogue-all "${CATALOGUE_ALL}" \
  --nsltp-evidence "${NSLTP_EVIDENCE}" \
  --out-dir "${PREP}"

QUERY=${PREP}/overlap_complete_start_stop_suborfs.faa
[[ -s "${QUERY}" ]] || { echo "No complete start-stop ORFs to audit" >&2; exit 1; }

"${DIAMOND}" blastp --query "${QUERY}" \
  --db "${SOURCE_RESULTS}/screen/references.dmnd" \
  --out "${OUT}/overlap_vs_labeled_AMP_refs.tsv" \
  --sensitive --evalue 0.001 --id 20 --max-target-seqs 0 --threads 2 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore

"${DIAMOND}" blastp --query "${QUERY}" \
  --db "${PROJECT}/intermediate/mcscanx/rso.dmnd" \
  --out "${OUT}/overlap_vs_rso_proteome.tsv" \
  --sensitive --evalue 0.001 --id 20 --max-target-seqs 25 --threads 2 \
  --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore

"${HMMSEARCH}" --cpu 2 --cut_ga --noali \
  --domtblout "${OUT}/overlap_family_domains.domtbl" \
  "${PROJECT}/data/hmm/AMP_models.hmm" "${QUERY}" \
  > "${OUT}/logs/overlap_hmmsearch.log"

"${BLASTP}" -task blastp-short \
  -query "${QUERY}" -subject "${NSLTP_LABELED_REFERENCE}" \
  -evalue 10 -seg no -comp_based_stats 0 \
  -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
  -out "${OUT}/overlap_vs_nsltp_labeled_reference.tsv"
python3 "${NSLTP_REFERENCE_SUMMARIZER}" \
  --blast "${OUT}/overlap_vs_nsltp_labeled_reference.tsv" \
  --queries "${QUERY}" \
  --output "${OUT}/overlap_nsltp_reference_similarity.tsv"

if [[ -f "${PREDGPI}/predgpi.py" ]]; then
  export PREDGPI_HOME="${PREDGPI}"
  python3 "${PREDGPI}/predgpi.py" -f "${QUERY}" \
    -o "${OUT}/overlap_start_stop_suborfs.predgpi.json" -m json \
    > "${OUT}/logs/overlap_predgpi.stdout.log" \
    2> "${OUT}/logs/overlap_predgpi.stderr.log"
else
  echo "PredGPI runtime unavailable" > "${OUT}/PREDGPI_NOT_RUN.txt"
fi

if [[ -x "${TMBED_PY}" && -d "${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc" ]]; then
  export HF_HOME="${REVISION}/08_reproducibility/models/huggingface_cache"
  export TRANSFORMERS_OFFLINE=1
  export HF_HUB_OFFLINE=1
  "${TMBED_PY}" -m tmbed predict \
    --fasta "${QUERY}" \
    --predictions "${OUT}/overlap_start_stop_suborfs.tmbed.pred" \
    --out-format 0 --no-use-gpu --threads 2 \
    --model-dir "${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc" \
    > "${OUT}/logs/overlap_tmbed.log" 2>&1
else
  echo "TMbed runtime unavailable" > "${OUT}/TMBED_NOT_RUN.txt"
fi

sha256sum "${SOURCE_RESULTS}/candidate_evidence.tsv" "${QUERY}" \
  "${CATALOGUE_ALL}" "${NSLTP_EVIDENCE}" "${NSLTP_LABELED_REFERENCE}" \
  > "${OUT}/input_identities.sha256"
date --iso-8601=seconds > "${OUT}/OVERLAP_FAMILY_AUDIT_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS overlap/altered family audit"
