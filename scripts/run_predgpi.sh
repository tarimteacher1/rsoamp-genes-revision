#!/usr/bin/env bash
set -euo pipefail

REV=/path/to/rsoamp/revision_R1_20260902
SOFT=$REV/08_reproducibility/software/predgpi
OUT=$REV/01_nslTP_curation
QUERY=$OUT/inputs/nsLTP_candidate_pool.faa
REF=$OUT/reference/positive_reviewed_arabidopsis_nsLTP.faa
mkdir -p "$OUT/predictions" "$OUT/logs"

export PREDGPI_HOME=$SOFT
python3 "$SOFT/predgpi.py" -f "$SOFT/testdata/test.fasta" \
    -o "$OUT/predictions/predgpi_fixture.json" -m json \
    > "$OUT/logs/predgpi_fixture.stdout.log" 2> "$OUT/logs/predgpi_fixture.stderr.log"
python3 "$SOFT/predgpi.py" -f "$QUERY" \
    -o "$OUT/predictions/nsLTP_candidates.predgpi.json" -m json \
    > "$OUT/logs/predgpi_candidates.stdout.log" 2> "$OUT/logs/predgpi_candidates.stderr.log"
python3 "$SOFT/predgpi.py" -f "$REF" \
    -o "$OUT/predictions/positive_reference.predgpi.json" -m json \
    > "$OUT/logs/predgpi_reference.stdout.log" 2> "$OUT/logs/predgpi_reference.stderr.log"

python3 - "$OUT" <<'PY'
import csv
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])

def parse(path):
    rows = []
    for record in json.load(path.open()):
        features = [feature for feature in record.get("features", []) if feature.get("description") == "GPI-anchor"]
        feature = features[0] if features else {}
        rows.append({
            "protein_id": record["accession"],
            "predgpi_GPI_anchor": "YES" if feature else "NO",
            "predgpi_omega_site": feature.get("begin", ""),
            "predgpi_score": feature.get("score", ""),
            "sequence_length": record["sequence"]["length"],
        })
    return rows

for source, destination in [
    ("nsLTP_candidates.predgpi.json", "nsLTP_candidates_predgpi.tsv"),
    ("positive_reference.predgpi.json", "positive_reference_predgpi.tsv"),
]:
    rows = parse(out / "predictions" / source)
    with (out / destination).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(destination, len(rows), "GPI", sum(row["predgpi_GPI_anchor"] == "YES" for row in rows))
PY

git -C "$SOFT" rev-parse HEAD > "$OUT/predgpi_git_commit.txt"
sha256sum "$SOFT/GPIDAT/"* "$QUERY" "$OUT/predictions/nsLTP_candidates.predgpi.json" \
    > "$OUT/predgpi_inputs_outputs.sha256"
