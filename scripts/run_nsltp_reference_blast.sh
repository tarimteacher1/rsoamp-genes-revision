#!/usr/bin/env bash
set -euo pipefail

REV=/path/to/rsoamp/revision_R1_20260902
OUT=$REV/01_nslTP_curation
REF=$OUT/reference/combined_labeled_reference.faa
QUERY=$OUT/inputs/nsLTP_candidate_pool.faa
DB=$OUT/reference/blastdb/combined_labeled_reference
BLAST=/path/to/user-home/tools/miniconda3/bin
mkdir -p "$OUT/reference/blastdb" "$OUT/predictions" "$OUT/logs"

"$BLAST/makeblastdb" -in "$REF" -dbtype prot -parse_seqids -out "$DB" \
    > "$OUT/logs/reference_makeblastdb.log" 2>&1
"$BLAST/blastp" -task blastp-short -query "$QUERY" -db "$DB" -evalue 10 \
    -seg no -comp_based_stats 0 -max_target_seqs 500 \
    -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
    -out "$OUT/predictions/nsLTP_candidates_vs_labeled_reference.tsv"

python3 - "$OUT" <<'PY'
import csv
import sys
from collections import defaultdict
from pathlib import Path

out = Path(sys.argv[1])
hits = defaultdict(list)
with (out / "predictions/nsLTP_candidates_vs_labeled_reference.tsv").open() as handle:
    for line in handle:
        q, s, pident, length, qlen, slen, evalue, bitscore = line.rstrip("\n").split("\t")
        parts = s.split("|")
        category = parts[1] if len(parts) >= 3 and parts[0] == "gnl" else parts[0]
        hits[q].append(
            {"subject": s, "category": category, "pident": float(pident), "alignment_length": int(length),
             "qlen": int(qlen), "slen": int(slen), "evalue": float(evalue), "bitscore": float(bitscore)}
        )

queries = []
current = None
for raw in (out / "inputs/nsLTP_candidate_pool.faa").open():
    if raw.startswith(">"):
        current = raw[1:].split()[0]
        queries.append(current)

rows = []
for query in queries:
    best = {}
    for category in ("POS", "NEG2S", "AMBPRO"):
        candidates = [hit for hit in hits.get(query, []) if hit["category"] == category]
        if candidates:
            best[category] = max(candidates, key=lambda hit: (hit["bitscore"], -hit["evalue"]))
    pos = best.get("POS", {})
    neg2s = best.get("NEG2S", {})
    ambpro = best.get("AMBPRO", {})
    alternatives = [hit for hit in (neg2s, ambpro) if hit]
    strongest_alternative = max(alternatives, key=lambda hit: hit["bitscore"]) if alternatives else {}
    margin = pos.get("bitscore", 0.0) - strongest_alternative.get("bitscore", 0.0)
    if pos and margin >= 10:
        interpretation = "positive_reference_preferred"
    elif neg2s and neg2s.get("bitscore", 0.0) - pos.get("bitscore", 0.0) >= 10 and neg2s.get("bitscore", 0.0) >= ambpro.get("bitscore", 0.0):
        interpretation = "explicit_2S_seed_storage_reference_preferred"
    elif ambpro and ambpro.get("bitscore", 0.0) - pos.get("bitscore", 0.0) >= 10:
        interpretation = "broad_prolamin_umbrella_reference_preferred"
    else:
        interpretation = "reference_similarity_ambiguous"
    rows.append({
        "protein_id": query,
        "best_positive_subject": pos.get("subject", ""),
        "best_positive_pident": pos.get("pident", ""),
        "best_positive_evalue": pos.get("evalue", ""),
        "best_positive_bitscore": pos.get("bitscore", ""),
        "best_explicit_2S_subject": neg2s.get("subject", ""),
        "best_explicit_2S_pident": neg2s.get("pident", ""),
        "best_explicit_2S_evalue": neg2s.get("evalue", ""),
        "best_explicit_2S_bitscore": neg2s.get("bitscore", ""),
        "best_ambiguous_prolamin_subject": ambpro.get("subject", ""),
        "best_ambiguous_prolamin_pident": ambpro.get("pident", ""),
        "best_ambiguous_prolamin_evalue": ambpro.get("evalue", ""),
        "best_ambiguous_prolamin_bitscore": ambpro.get("bitscore", ""),
        "positive_minus_strongest_alternative_bitscore": f"{margin:.1f}",
        "reference_similarity_interpretation": interpretation,
    })

path = out / "nsLTP_labeled_reference_similarity.tsv"
with path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print("Wrote", len(rows), "rows to", path)
PY

sha256sum "$REF" "$QUERY" "$OUT/predictions/nsLTP_candidates_vs_labeled_reference.tsv" \
    > "$OUT/reference_blast_outputs.sha256"
