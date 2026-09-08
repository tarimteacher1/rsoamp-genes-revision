#!/usr/bin/env python3
"""Summarize short-protein BLAST hits against labeled nsLTP/prolamin references."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def fasta_ids(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]


def main() -> None:
    args = parse_args()
    hits = defaultdict(list)
    with args.blast.open(encoding="utf-8") as handle:
        for line in handle:
            query, subject, pident, length, qlen, slen, evalue, bitscore = line.rstrip("\n").split("\t")
            parts = subject.split("|")
            category = parts[1] if len(parts) >= 3 and parts[0] == "gnl" else parts[0]
            hits[query].append(
                {
                    "subject": subject, "category": category, "pident": float(pident),
                    "alignment_length": int(length), "qlen": int(qlen), "slen": int(slen),
                    "evalue": float(evalue), "bitscore": float(bitscore),
                }
            )

    rows = []
    for query in fasta_ids(args.queries):
        best = {}
        for category in ("POS", "NEG2S", "AMBPRO"):
            candidates = [hit for hit in hits.get(query, []) if hit["category"] == category]
            if candidates:
                best[category] = max(candidates, key=lambda hit: (hit["bitscore"], -hit["evalue"]))
        pos, neg2s, ambpro = best.get("POS", {}), best.get("NEG2S", {}), best.get("AMBPRO", {})
        alternatives = [hit for hit in (neg2s, ambpro) if hit]
        strongest = max(alternatives, key=lambda hit: hit["bitscore"]) if alternatives else {}
        margin = pos.get("bitscore", 0.0) - strongest.get("bitscore", 0.0)
        if pos and margin >= 10:
            interpretation = "positive_reference_preferred"
        elif neg2s and neg2s.get("bitscore", 0.0) - pos.get("bitscore", 0.0) >= 10 and neg2s.get("bitscore", 0.0) >= ambpro.get("bitscore", 0.0):
            interpretation = "explicit_2S_seed_storage_reference_preferred"
        elif ambpro and ambpro.get("bitscore", 0.0) - pos.get("bitscore", 0.0) >= 10:
            interpretation = "broad_prolamin_umbrella_reference_preferred"
        else:
            interpretation = "reference_similarity_ambiguous"
        rows.append(
            {
                "protein_id": query,
                "best_positive_subject": pos.get("subject", ""),
                "best_positive_pident": pos.get("pident", ""),
                "best_positive_evalue": pos.get("evalue", ""),
                "best_positive_bitscore": pos.get("bitscore", ""),
                "best_explicit_2S_subject": neg2s.get("subject", ""),
                "best_explicit_2S_bitscore": neg2s.get("bitscore", ""),
                "best_ambiguous_prolamin_subject": ambpro.get("subject", ""),
                "best_ambiguous_prolamin_bitscore": ambpro.get("bitscore", ""),
                "positive_minus_strongest_alternative_bitscore": f"{margin:.1f}",
                "reference_similarity_interpretation": interpretation,
            }
        )
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} labeled-reference summaries")


if __name__ == "__main__":
    main()
