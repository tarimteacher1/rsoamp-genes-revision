#!/usr/bin/env python3
"""Audit R1.M6 priority models against frozen annotation and old rescue loci."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Refusing to write an empty neighborhood audit")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_gtf(path: Path) -> dict[str, dict[str, object]]:
    transcripts: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = dict(ATTRIBUTE_RE.findall(fields[8]))
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            start, end = int(fields[3]), int(fields[4])
            model = transcripts.setdefault(
                tid,
                {
                    "chrom": fields[0],
                    "strand": fields[6],
                    "start": start,
                    "end": end,
                    "gene_id": attrs.get("gene_id", ""),
                    "exons": [],
                },
            )
            model["start"] = min(int(model["start"]), start)
            model["end"] = max(int(model["end"]), end)
            if fields[2] == "exon":
                model["exons"].append((start, end))
    for model in transcripts.values():
        model["exons"] = sorted(set(model["exons"]))
    return transcripts


def interval_distance(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    if start_a <= end_b and start_b <= end_a:
        return 0
    if end_a < start_b:
        return start_b - end_a - 1
    return start_a - end_b - 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-audit", required=True, type=Path)
    parser.add_argument("--candidate-gtf", required=True, type=Path)
    parser.add_argument("--reference-gtf", required=True, type=Path)
    parser.add_argument("--old-rescue", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    candidates = read_tsv(args.candidate_audit)
    candidate_models = parse_gtf(args.candidate_gtf)
    reference_models = parse_gtf(args.reference_gtf)
    old_rescue = read_tsv(args.old_rescue)

    reference_by_chrom: dict[str, list[tuple[str, dict[str, object]]]] = defaultdict(list)
    for tid, model in reference_models.items():
        reference_by_chrom[str(model["chrom"])].append((tid, model))

    old_by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in old_rescue:
        old_by_chrom[row["chromosome"]].append(row)

    exon_owners: dict[tuple[str, str, int, int], set[str]] = defaultdict(set)
    for tid, model in candidate_models.items():
        for start, end in model["exons"]:
            exon_owners[(str(model["chrom"]), str(model["strand"]), start, end)].add(tid)

    peptide_owners: dict[str, list[str]] = defaultdict(list)
    for row in candidates:
        peptide_owners[row["protein_sequence"]].append(row["transcript_id"])

    output: list[dict[str, object]] = []
    for row in candidates:
        tid = row["transcript_id"]
        model = candidate_models[tid]
        chrom = str(model["chrom"])
        start, end = int(model["start"]), int(model["end"])
        exons = list(model["exons"])
        introns = [exons[index + 1][0] - exons[index][1] - 1 for index in range(len(exons) - 1)]

        reference_distances = [
            (interval_distance(start, end, int(other["start"]), int(other["end"])), other_tid, other)
            for other_tid, other in reference_by_chrom.get(chrom, [])
        ]
        nearest_distance, nearest_tid, nearest = min(reference_distances, default=(-1, "", {}), key=lambda item: item[0])
        overlapping_reference = sorted(
            other_tid for distance, other_tid, _ in reference_distances if distance == 0
        )

        rescue_distances = [
            (
                interval_distance(
                    start,
                    end,
                    int(hit["cluster_start_1based"]),
                    int(hit["cluster_end_1based"]),
                ),
                hit,
            )
            for hit in old_by_chrom.get(chrom, [])
        ]
        old_distance, nearest_old = min(rescue_distances, default=(-1, {}), key=lambda item: item[0])
        overlapping_old = sorted(hit["screen_locus_id"] for distance, hit in rescue_distances if distance == 0)

        shared_exon_models = sorted(
            {
                other_tid
                for exon_start, exon_end in exons
                for other_tid in exon_owners[(chrom, str(model["strand"]), exon_start, exon_end)]
                if other_tid != tid
            }
        )
        duplicate_peptide_models = sorted(other for other in peptide_owners[row["protein_sequence"]] if other != tid)
        structure_flags = []
        if introns and max(introns) > 10000:
            structure_flags.append("INTRON_GT_10KB")
        if shared_exon_models:
            structure_flags.append("EXACT_EXON_REUSE")
        if duplicate_peptide_models:
            structure_flags.append("IDENTICAL_PROTEIN_AT_OTHER_MODEL")

        output.append(
            {
                "provisional_locus_id": row["provisional_locus_id"],
                "transcript_id": tid,
                "chrom": chrom,
                "model_start": start,
                "model_end": end,
                "strand": model["strand"],
                "model_span_bp": end - start + 1,
                "exon_count": len(exons),
                "intron_lengths_bp": ";".join(map(str, introns)) if introns else "NA",
                "maximum_intron_bp": max(introns) if introns else 0,
                "overlapping_frozen_transcript_count": len(overlapping_reference),
                "overlapping_frozen_transcript_ids": ";".join(overlapping_reference),
                "nearest_frozen_transcript_id": nearest_tid,
                "nearest_frozen_gene_id": nearest.get("gene_id", "") if nearest else "",
                "nearest_frozen_transcript_distance_bp": nearest_distance,
                "nearest_frozen_transcript_coordinates": (
                    f"{nearest.get('chrom')}:{nearest.get('start')}-{nearest.get('end')}({nearest.get('strand')})"
                    if nearest
                    else ""
                ),
                "overlapping_old_rescue_locus_count": len(overlapping_old),
                "overlapping_old_rescue_locus_ids": ";".join(overlapping_old),
                "nearest_old_rescue_locus_id": nearest_old.get("screen_locus_id", "") if nearest_old else "",
                "nearest_old_rescue_distance_bp": old_distance,
                "candidate_models_sharing_exact_exons": ";".join(shared_exon_models),
                "candidate_models_with_identical_protein": ";".join(duplicate_peptide_models),
                "candidate_structure_review_flags": ";".join(structure_flags) if structure_flags else "NONE",
                "interpretation": "Interval comparison only; transcript and family acceptance require independent read/model evidence",
            }
        )

    write_tsv(args.output, output)
    print(f"Wrote {len(output)} candidate neighborhood audit rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
