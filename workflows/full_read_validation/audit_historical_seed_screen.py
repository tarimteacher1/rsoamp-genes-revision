#!/usr/bin/env python3
"""Reconcile R1.M6 candidates with the frozen whole-genome seed screen."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_seed_counts(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with path.open(encoding="ascii") as handle:
        for raw in handle:
            if raw.startswith(">"):
                counts[raw[1:].split("|", 1)[0].strip()] += 1
    return counts


def read_historical_hits(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise ValueError(f"Unexpected historical tblastn row: {raw.rstrip()}")
            query, chrom, pident, length, evalue, bitscore, sstart, send = fields
            rows.append(
                {
                    "query": query,
                    "chrom": chrom,
                    "pident": float(pident),
                    "length": int(length),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                    "start": min(int(sstart), int(send)),
                    "end": max(int(sstart), int(send)),
                }
            )
    return rows


def read_candidate_seed_hits(path: Path) -> dict[str, list[dict[str, object]]]:
    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise ValueError(f"Unexpected candidate/seed BLAST row: {raw.rstrip()}")
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore = fields
            hits[qseqid.removesuffix("|start_stop")].append(
                {
                    "subject": sseqid,
                    "pident": float(pident),
                    "length": int(length),
                    "qlen": int(qlen),
                    "slen": int(slen),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                }
            )
    for records in hits.values():
        records.sort(key=lambda row: (-float(row["bitscore"]), float(row["evalue"])))
    return hits


def infer_family(row: dict[str, str]) -> str:
    evidence = f'{row.get("family_HMMs", "")} {row.get("best_AMP_reference", "")}'.lower()
    if "pf00304" in evidence or "defensin" in evidence or "rsdef" in evidence:
        return "Defensin"
    if "pf00234" in evidence or "nsltp" in evidence or "rsltp" in evidence:
        return "nsLTP"
    return "Unresolved"


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-audit", required=True, type=Path)
    parser.add_argument("--sequence-manifest", required=True, type=Path)
    parser.add_argument("--historical-raw-hits", required=True, type=Path)
    parser.add_argument("--candidate-vs-historical-seeds", required=True, type=Path)
    parser.add_argument("--historical-seeds", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    candidates = read_tsv(args.candidate_audit)
    manifest = {row["transcript_id"]: row for row in read_tsv(args.sequence_manifest)}
    historical = read_historical_hits(args.historical_raw_hits)
    seed_hits = read_candidate_seed_hits(args.candidate_vs_historical_seeds)
    seed_counts = read_seed_counts(args.historical_seeds)

    if len(candidates) != 10 or set(manifest) & {row["transcript_id"] for row in candidates} != {
        row["transcript_id"] for row in candidates
    }:
        raise ValueError("Priority candidate identities do not join 10/10")
    if sum(seed_counts.values()) != 128:
        raise ValueError(f"Expected 128 historical seeds, observed {sum(seed_counts.values())}")

    rows: list[dict[str, object]] = []
    for candidate in candidates:
        tid = candidate["transcript_id"]
        sequence = manifest[tid]
        family = infer_family(sequence)
        start = int(candidate["model_start"])
        end = int(candidate["model_end"])
        overlapping = [
            hit
            for hit in historical
            if hit["chrom"] == candidate["chrom"]
            and overlaps(start, end, int(hit["start"]), int(hit["end"]))
        ]
        best = seed_hits.get(candidate["orf_id"], [{}])[0]
        best_evalue = best.get("evalue")
        if family == "nsLTP" and seed_counts.get("nsLTP", 0) == 0:
            interpretation = "HISTORICAL_PANEL_CONTAINED_NO_NSLTP_SEEDS"
        elif not overlapping and best_evalue is not None and float(best_evalue) > 1e-5:
            interpretation = "NO_FROZEN_HSP;DIVERGENT_FROM_HISTORICAL_SEEDS_AT_OLD_THRESHOLD"
        elif not overlapping:
            interpretation = "NO_FROZEN_HSP;SEARCH_SPACE_OR_SCORING_DIFFERENCE_REQUIRES_CAUTION"
        else:
            interpretation = "OVERLAPS_FROZEN_HISTORICAL_HSP"
        rows.append(
            {
                "transcript_id": tid,
                "orf_id": candidate["orf_id"],
                "family_hint": family,
                "chrom": candidate["chrom"],
                "model_start": start,
                "model_end": end,
                "historical_seed_count_for_family": seed_counts.get(family, 0),
                "historical_raw_HSP_overlap_count": len(overlapping),
                "best_historical_seed": best.get("subject", ""),
                "best_candidate_to_seed_pident": best.get("pident", "NA"),
                "best_candidate_to_seed_alignment_aa": best.get("length", "NA"),
                "best_candidate_to_seed_query_length_aa": best.get("qlen", "NA"),
                "best_candidate_to_seed_evalue_blastp_short": best_evalue if best_evalue is not None else "NA",
                "best_candidate_to_seed_bitscore_blastp_short": best.get("bitscore", "NA"),
                "frozen_tblastn_evalue_threshold": "1e-5",
                "interpretation": interpretation,
                "comparability_note": (
                    "Candidate-to-seed blastp-short and seed-to-whole-genome tblastn E-values are not "
                    "numerically interchangeable; this audit explains panel coverage and observed non-overlap."
                ),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "historical_seed_total": sum(seed_counts.values()),
        "historical_seed_family_counts": dict(sorted(seed_counts.items())),
        "priority_candidates": len(rows),
        "priority_candidates_overlapping_historical_raw_HSPs": sum(
            int(row["historical_raw_HSP_overlap_count"]) > 0 for row in rows
        ),
        "interpretation": (
            "The frozen 128-seed screen was not exhaustive for the new R1.M6 candidates. Its result must be "
            "reported as conditional on the historical panel and threshold, not as proof that no other AMP-like "
            "locus exists."
        ),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
