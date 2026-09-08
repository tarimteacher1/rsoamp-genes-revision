#!/usr/bin/env bash
set -euo pipefail

REVISION=${REVISION:-/path/to/rsoamp/revision_R1_20260902}
OUT=${REVISION}/05_comparative_genomics/mcscanx_rso_tamarix
DIAMOND=${DIAMOND:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/diamond}
MCSCANX=${MCSCANX:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/MCScanX}
SAMTOOLS=${SAMTOOLS:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/samtools}
THREADS=${THREADS:-48}
mkdir -p "${OUT}" "${REVISION}/05_comparative_genomics/logs" "${REVISION}/08_reproducibility/commands"

TAU_GENOME=${REVISION}/05_comparative_genomics/inputs/Tamarix_austromongolica/tau_genome.fasta
if [[ ! -s "${TAU_GENOME}.fai" ]]; then
  "${SAMTOOLS}" faidx "${TAU_GENOME}"
fi

python3 "${REVISION}/05_comparative_genomics/prepare_rso_tamarix_mcscanx.py" prepare
cd "${OUT}"

"${DIAMOND}" makedb --in rso_tau.faa -d rso_tau \
  > "${REVISION}/05_comparative_genomics/logs/rso_tau_diamond_makedb.log" 2>&1
"${DIAMOND}" blastp -q rso_tau.faa -d rso_tau -o rso_tau.blast \
  --outfmt 6 --evalue 1e-5 --max-target-seqs 10 --threads "${THREADS}" --sensitive \
  > "${REVISION}/05_comparative_genomics/logs/rso_tau_diamond_blastp.stdout.log" \
  2> "${REVISION}/05_comparative_genomics/logs/rso_tau_diamond_blastp.stderr.log"

"${MCSCANX}" rso_tau \
  > "${REVISION}/05_comparative_genomics/logs/rso_tau_mcscanx.stdout.log" \
  2> "${REVISION}/05_comparative_genomics/logs/rso_tau_mcscanx.stderr.log"
[[ -s rso_tau.collinearity ]] || { echo "MCScanX collinearity output is empty" >&2; exit 1; }

python3 "${REVISION}/05_comparative_genomics/prepare_rso_tamarix_mcscanx.py" parse
sha256sum rso_tau.faa rso_tau.gff rso_tau.blast rso_tau.collinearity \
  > "${REVISION}/05_comparative_genomics/rso_tamarix_mcscanx.sha256"
date --iso-8601=seconds > "${REVISION}/05_comparative_genomics/RSO_TAMARIX_MCSCANX_COMPLETE.PASS"
echo "PASS R. soongarica-T. austromongolica MCScanX comparison"
