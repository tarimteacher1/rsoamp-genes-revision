#!/usr/bin/env bash
set -euo pipefail

ORIGINAL=${ORIGINAL:-/path/to/rsoamp}
REVISION=${REVISION:-${ORIGINAL}/revision_R1_20260902}
MAKEBLASTDB=${MAKEBLASTDB:-/path/to/user-home/tools/miniconda3/bin/makeblastdb}
TBLASTN=${TBLASTN:-/path/to/user-home/tools/miniconda3/bin/tblastn}
MINIPROT=${MINIPROT:-/path/to/user-home/tools/miniconda3/envs/busco1/bin/miniprot}
THREADS=${THREADS:-32}

OUT=${REVISION}/05_comparative_genomics
IN=${OUT}/inputs/Tamarix_austromongolica
PRED=${OUT}/predictions
LOG=${OUT}/logs
DB=${OUT}/blastdb/tau_genome
GENOME=${IN}/tau_genome.fasta
QUERY=${PRED}/rso_annotated_defensins.faa
mkdir -p "$(dirname "${DB}")" "${PRED}" "${LOG}" "${REVISION}/08_reproducibility/commands"

echo '566c148566390c66c677207303f3b305  '"${GENOME}" | md5sum -c -
python3 "${OUT}/prepare_tamarix_defensin_rescue.py" prepare

if [[ ! -s "${DB}.nsq" ]]; then
  "${MAKEBLASTDB}" -in "${GENOME}" -dbtype nucl -parse_seqids -out "${DB}" \
    > "${LOG}/tamarix_makeblastdb.log" 2>&1
fi

cat > "${REVISION}/08_reproducibility/commands/tamarix_defensin_tblastn.command.txt" <<EOF
${TBLASTN} -query ${QUERY} -db ${DB} -evalue 1e-4 -seg no -max_hsps 50 -max_target_seqs 500 -num_threads ${THREADS} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'
EOF
"${TBLASTN}" -query "${QUERY}" -db "${DB}" -evalue 1e-4 -seg no \
  -max_hsps 50 -max_target_seqs 500 -num_threads "${THREADS}" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' \
  -out "${PRED}/tamarix_defensin_tblastn.tsv" \
  > "${LOG}/tamarix_defensin_tblastn.stdout.log" \
  2> "${LOG}/tamarix_defensin_tblastn.stderr.log"

"${MINIPROT}" -t "${THREADS}" --gff "${GENOME}" "${QUERY}" \
  > "${PRED}/tamarix_defensin_miniprot.gff" \
  2> "${LOG}/tamarix_defensin_miniprot.log" || true
"${MINIPROT}" -t "${THREADS}" --trans "${GENOME}" "${QUERY}" \
  > "${PRED}/tamarix_defensin_miniprot_translated.txt" \
  2> "${LOG}/tamarix_defensin_miniprot_translated.log" || true

python3 "${OUT}/prepare_tamarix_defensin_rescue.py" parse
sha256sum "${GENOME}" "${QUERY}" "${PRED}/tamarix_defensin_tblastn.tsv" \
  "${PRED}/tamarix_defensin_miniprot.gff" "${PRED}/tamarix_defensin_miniprot_translated.txt" \
  > "${OUT}/tamarix_defensin_rescue.sha256"
date --iso-8601=seconds > "${OUT}/TAMARIX_DEFENSIN_RESCUE_COMPLETE.PASS"
echo "PASS Tamarix defensin genome rescue"
