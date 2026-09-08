#!/usr/bin/env bash
set -euo pipefail

ORIGINAL=${ORIGINAL:-/path/to/rsoamp}
REVISION=${REVISION:-${ORIGINAL}/revision_R1_20260902}
PYTHON=${PYTHON:-/path/to/user-home/tools/miniconda3/envs/atac_proc/bin/python}
THREADS=${THREADS:-16}

SOFTWARE=${REVISION}/08_reproducibility/software/TMbed
ENV=${REVISION}/08_reproducibility/envs/tmbed
MODEL=${REVISION}/08_reproducibility/models/prot_t5_xl_half_uniref50_enc
INPUT=${REVISION}/01_nslTP_curation/inputs/nsLTP_candidate_pool.faa
OUTPUT=${REVISION}/01_nslTP_curation/predictions/nsLTP_candidates.tmbed.pred
SUMMARY=${REVISION}/01_nslTP_curation/nsLTP_candidates_tmbed.tsv
LOG=${REVISION}/01_nslTP_curation/logs/tmbed.log
mkdir -p "$(dirname "${OUTPUT}")" "$(dirname "${LOG}")" "${MODEL}" "$(dirname "${ENV}")"

exec > >(tee "${LOG}") 2>&1
echo "[$(date --iso-8601=seconds)] START TMbed"
[[ -s "${INPUT}" ]] || { echo "Missing input: ${INPUT}" >&2; exit 1; }

if [[ ! -x "${ENV}/bin/python" ]]; then
  "${PYTHON}" -m venv "${ENV}"
fi
"${ENV}/bin/python" -m pip install --upgrade 'pip<26'
"${ENV}/bin/python" -m pip install \
  --index-url https://download.pytorch.org/whl/cpu \
  'torch==2.5.1'
"${ENV}/bin/python" -m pip install \
  'h5py==3.12.1' \
  'numpy==1.26.4' \
  'sentencepiece==0.2.0' \
  'tqdm==4.67.1' \
  'transformers==4.46.3' \
  'typer==0.4.2' \
  'click==8.1.7'
"${ENV}/bin/python" -m pip install --no-deps "${SOFTWARE}"

export HF_HOME=${REVISION}/08_reproducibility/models/huggingface_cache
export HF_ENDPOINT=${HF_ENDPOINT:-https://hf-mirror.com}
"${ENV}/bin/python" -m tmbed predict \
  --fasta "${INPUT}" \
  --predictions "${OUTPUT}" \
  --out-format 0 \
  --no-use-gpu \
  --threads "${THREADS}" \
  --model-dir "${MODEL}"

"${ENV}/bin/python" "${REVISION}/01_nslTP_curation/parse_tmbed.py" "${OUTPUT}" "${SUMMARY}"
[[ "$(awk 'END {print NR-1}' "${SUMMARY}")" -eq 54 ]] || {
  echo "Expected 54 TMbed rows" >&2
  exit 1
}

"${ENV}/bin/python" -m pip freeze > "${REVISION}/08_reproducibility/tmbed_environment_freeze.txt"
git -C "${SOFTWARE}" rev-parse HEAD > "${REVISION}/08_reproducibility/tmbed_git_commit.txt"
date --iso-8601=seconds > "${REVISION}/01_nslTP_curation/TMBED_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS TMbed"
