#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/rsoamp/revision_R1_20260902}"
DIR="$BASE/05_comparative_genomics"
OUT="$DIR/tamarix_nsltp_recuration"
TAU="$DIR/inputs/Tamarix_austromongolica/tau.longest_pep.fasta"
PRELIM="$DIR/tamarix_amp_evidence_preliminary.tsv"
REF="$BASE/01_nslTP_curation/reference/combined_labeled_reference.faa"
BLAST_BIN="${BLAST_BIN:-/path/to/user-home/tools/miniconda3/bin}"
PREDGPI="$BASE/08_reproducibility/software/predgpi"
TMBED_ENV="$BASE/08_reproducibility/envs/tmbed"
TMBED_MODEL="$BASE/08_reproducibility/models/prot_t5_xl_half_uniref50_enc"
mkdir -p "$OUT/inputs" "$OUT/predictions" "$OUT/logs" "$OUT/reference"

export HF_HOME="$BASE/08_reproducibility/models/huggingface_cache"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

python3 - "$PRELIM" "$TAU" "$OUT/inputs/nsLTP_candidate_pool.faa" <<'PY'
import csv, sys
from pathlib import Path
prelim, fasta, output = map(Path, sys.argv[1:])
with prelim.open() as handle:
    ids = {row['protein_id'] for row in csv.DictReader(handle, delimiter='\t') if 'nsLTP' in row['hmm_families'].split(';')}
records, current, chunks = {}, None, []
for line in fasta.read_text().splitlines():
    if line.startswith('>'):
        if current is not None: records[current] = ''.join(chunks)
        current, chunks = line[1:].split()[0], []
    else: chunks.append(line.strip())
if current is not None: records[current] = ''.join(chunks)
if ids - records.keys(): raise SystemExit('Missing Tamarix candidate sequences')
with output.open('w') as handle:
    for protein_id in sorted(ids): handle.write(f'>{protein_id}\n{records[protein_id]}\n')
print('candidate_pool', len(ids))
PY

QUERY="$OUT/inputs/nsLTP_candidate_pool.faa"
DB="$OUT/reference/combined_labeled_reference"
"$BLAST_BIN/makeblastdb" -in "$REF" -dbtype prot -parse_seqids -out "$DB" > "$OUT/logs/makeblastdb.log" 2>&1
"$BLAST_BIN/blastp" -task blastp-short -query "$QUERY" -db "$DB" -evalue 10 -seg no -comp_based_stats 0 \
    -max_target_seqs 500 -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
    -out "$OUT/predictions/candidates_vs_labeled_reference.tsv"
python3 "$DIR/summarize_labeled_reference_blast.py" \
    --blast "$OUT/predictions/candidates_vs_labeled_reference.tsv" --queries "$QUERY" \
    --output "$OUT/nsLTP_labeled_reference_similarity.tsv"

export PREDGPI_HOME="$PREDGPI"
python3 "$PREDGPI/predgpi.py" -f "$QUERY" -o "$OUT/predictions/candidates.predgpi.json" -m json \
    > "$OUT/logs/predgpi.stdout.log" 2> "$OUT/logs/predgpi.stderr.log"
python3 - "$OUT/predictions/candidates.predgpi.json" "$OUT/nsLTP_candidates_predgpi.tsv" <<'PY'
import csv, json, sys
rows=[]
for record in json.load(open(sys.argv[1])):
    features=[f for f in record.get('features',[]) if f.get('description')=='GPI-anchor']
    feature=features[0] if features else {}
    rows.append({'protein_id':record['accession'],'predgpi_GPI_anchor':'YES' if feature else 'NO','predgpi_omega_site':feature.get('begin',''),'predgpi_score':feature.get('score','')})
with open(sys.argv[2],'w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
PY

"$TMBED_ENV/bin/python" -m tmbed predict --fasta "$QUERY" \
    --predictions "$OUT/predictions/candidates.tmbed.pred" --out-format 0 --no-use-gpu --threads 16 --model-dir "$TMBED_MODEL" \
    > "$OUT/logs/tmbed.stdout.log" 2> "$OUT/logs/tmbed.stderr.log"
"$TMBED_ENV/bin/python" "$BASE/01_nslTP_curation/parse_tmbed.py" \
    "$OUT/predictions/candidates.tmbed.pred" "$OUT/nsLTP_candidates_tmbed.tsv"

python3 "$DIR/build_tamarix_nsltp_final.py" \
    --preliminary "$PRELIM" --fasta "$QUERY" \
    --predgpi "$OUT/nsLTP_candidates_predgpi.tsv" --tmbed "$OUT/nsLTP_candidates_tmbed.tsv" \
    --similarity "$OUT/nsLTP_labeled_reference_similarity.tsv" \
    --output "$DIR/tamarix_nsLTP_evidence_matrix.tsv" \
    --canonical-fasta "$DIR/predictions/tamarix_nsLTP_final_canonical.faa"

sha256sum "$QUERY" "$REF" "$DIR/tamarix_nsLTP_evidence_matrix.tsv" \
    > "$DIR/tamarix_nsLTP_recuration.sha256"
printf 'PASS\tTamarix nsLTP candidates processed with the same non-phylogenetic decision gates\n' \
    > "$DIR/TAMARIX_NSLTP_RECURATION_COMPLETE.PASS"
