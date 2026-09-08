#!/usr/bin/env python3
"""Audit fixed-threshold recursive genome hits from R1.M6 candidate proteins."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_hits(path: Path) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 14:
                raise ValueError(f"Unexpected tblastn row: {raw.rstrip()}")
            qseqid, chrom, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = fields
            left, right = sorted((int(sstart), int(send)))
            hits.append(
                {
                    "query": qseqid.removesuffix("|start_stop"),
                    "chrom": chrom,
                    "strand": "+" if int(sstart) <= int(send) else "-",
                    "start": left,
                    "end": right,
                    "qstart": min(int(qstart), int(qend)),
                    "qend": max(int(qstart), int(qend)),
                    "qlen": int(qlen),
                    "pident": float(pident),
                    "length": int(length),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                }
            )
    return hits


def parse_gtf(path: Path) -> list[dict[str, object]]:
    models: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = dict(ATTRIBUTE_RE.findall(fields[8]))
            tid = attrs.get("transcript_id") or attrs.get("ID")
            if not tid:
                continue
            model = models.setdefault(
                tid,
                {
                    "id": tid,
                    "gene_id": attrs.get("gene_id", ""),
                    "chrom": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                },
            )
            model["start"] = min(int(model["start"]), int(fields[3]))
            model["end"] = max(int(model["end"]), int(fields[4]))
    return list(models.values())


def merge_intervals(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if overlaps(a_start, a_end, b_start, b_end):
        return 0
    return b_start - a_end if a_end < b_start else a_start - b_end


def cluster_hits(hits: list[dict[str, object]], max_gap: int) -> list[list[dict[str, object]]]:
    groups: list[list[dict[str, object]]] = []
    for key in sorted({(str(hit["chrom"]), str(hit["strand"])) for hit in hits}):
        selected = sorted(
            (hit for hit in hits if (str(hit["chrom"]), str(hit["strand"])) == key),
            key=lambda hit: (int(hit["start"]), int(hit["end"])),
        )
        current: list[dict[str, object]] = []
        current_end = -1
        for hit in selected:
            if current and int(hit["start"]) > current_end + max_gap:
                groups.append(current)
                current = []
                current_end = -1
            current.append(hit)
            current_end = max(current_end, int(hit["end"])) if current_end >= 0 else int(hit["end"])
        if current:
            groups.append(current)
    return groups


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--priority-models", required=True, type=Path)
    parser.add_argument("--candidate-gtf", required=True, type=Path)
    parser.add_argument("--reference-gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--window-bed", type=Path)
    parser.add_argument("--window-extension", type=int, default=2000)
    parser.add_argument("--cluster-gap", type=int, default=2000)
    args = parser.parse_args()

    hits = parse_hits(args.hits)
    priority = read_tsv(args.priority_models)
    candidates = parse_gtf(args.candidate_gtf)
    reference = parse_gtf(args.reference_gtf)
    clusters = cluster_hits(hits, args.cluster_gap)
    rows: list[dict[str, object]] = []

    for number, cluster in enumerate(clusters, 1):
        chrom = str(cluster[0]["chrom"])
        strand = str(cluster[0]["strand"])
        start = min(int(hit["start"]) for hit in cluster)
        end = max(int(hit["end"]) for hit in cluster)
        query_hits: dict[str, list[dict[str, object]]] = defaultdict(list)
        for hit in cluster:
            query_hits[str(hit["query"])].append(hit)
        query_stats = []
        for query, records in query_hits.items():
            qlen = int(records[0]["qlen"])
            covered = merge_intervals([(int(record["qstart"]), int(record["qend"])) for record in records])
            query_stats.append(
                (
                    query,
                    covered / qlen,
                    sum(float(record["bitscore"]) for record in records),
                    max(float(record["pident"]) for record in records),
                )
            )
        query_stats.sort(key=lambda item: (-item[1], -item[2], item[0]))
        best_query, best_coverage, best_bitscore_sum, best_identity = query_stats[0]
        priority_overlaps = sorted(
            row["transcript_id"]
            for row in priority
            if row["chrom"] == chrom
            and row["strand"] == strand
            and overlaps(start, end, int(row["model_start"]), int(row["model_end"]))
        )
        candidate_overlaps = sorted(
            str(model["id"])
            for model in candidates
            if model["chrom"] == chrom
            and model["strand"] == strand
            and overlaps(start, end, int(model["start"]), int(model["end"]))
        )
        reference_overlaps = sorted(
            str(model["id"])
            for model in reference
            if model["chrom"] == chrom
            and model["strand"] == strand
            and overlaps(start, end, int(model["start"]), int(model["end"]))
        )
        nearest_candidate_distance = min(
            (
                distance(start, end, int(model["start"]), int(model["end"]))
                for model in candidates
                if model["chrom"] == chrom and model["strand"] == strand
            ),
            default=-1,
        )
        if priority_overlaps:
            interpretation = "EXPECTED_PRIORITY_CANDIDATE_LOCUS"
        elif candidate_overlaps:
            interpretation = "OTHER_DISCOVERY_MODEL_LOCUS"
        elif reference_overlaps:
            interpretation = "FROZEN_ANNOTATED_LOCUS"
        else:
            interpretation = "UNANNOTATED_GENOME_HOMOLOGY_CLUSTER_REQUIRES_MODEL_AND_READ_AUDIT"
        rows.append(
            {
                "recursive_cluster_id": f"R1M6_REC{number:03d}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "span_bp": end - start + 1,
                "HSP_count": len(cluster),
                "query_count": len(query_hits),
                "queries": ";".join(sorted(query_hits)),
                "best_query": best_query,
                "best_query_aggregate_coverage_fraction": f"{best_coverage:.6f}",
                "best_query_bitscore_sum": f"{best_bitscore_sum:.1f}",
                "best_HSP_identity_percent": f"{best_identity:.3f}",
                "minimum_evalue": f"{min(float(hit['evalue']) for hit in cluster):.3g}",
                "overlapping_priority_models": ";".join(priority_overlaps),
                "overlapping_discovery_models": ";".join(candidate_overlaps),
                "overlapping_frozen_transcripts": ";".join(reference_overlaps),
                "nearest_same_strand_discovery_model_distance_bp": nearest_candidate_distance,
                "interpretation": interpretation,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[str(row["interpretation"])] += 1
    summary = {
        "query_HSPs": len(hits),
        "recursive_clusters": len(rows),
        "cluster_gap_bp": args.cluster_gap,
        "interpretation_counts": dict(sorted(counts.items())),
        "unannotated_clusters": [
            row["recursive_cluster_id"]
            for row in rows
            if str(row["interpretation"]).startswith("UNANNOTATED_")
        ],
        "scope": "Fixed candidate-seed tblastn E < 1e-5 screen; a homology cluster is not a gene call.",
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if args.window_bed:
        with args.window_bed.open("w", encoding="utf-8", newline="") as handle:
            for row in rows:
                if not str(row["interpretation"]).startswith("UNANNOTATED_"):
                    continue
                start0 = max(0, int(row["start"]) - 1 - args.window_extension)
                end0 = int(row["end"]) + args.window_extension
                handle.write(
                    f'{row["chrom"]}\t{start0}\t{end0}\t{row["recursive_cluster_id"]}|{row["strand"]}\n'
                )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
