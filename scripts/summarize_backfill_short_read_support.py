#!/usr/bin/env python3
"""Summarize exact-locus short-read support from target-filtered STAR BAMs."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--bam-root", required=True, type=Path)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--per-sample-output", required=True, type=Path)
    parser.add_argument("--locus-output", required=True, type=Path)
    return parser.parse_args()


def parse_loci(path: Path) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 6:
        raise ValueError(f"Expected six backfill loci, found {len(rows)}")
    parsed = []
    for row in rows:
        chrom = row["historical_locus"].split(":", 1)[0]
        start = min(int(row["coding_start"]), int(row["coding_end"]))
        end = max(int(row["coding_start"]), int(row["coding_end"]))
        parsed.append({**row, "chrom": chrom, "coding_lo": start, "coding_hi": end})
    return parsed


def run_text(command: list[str]) -> str:
    return subprocess.run(command, check=True, text=True, capture_output=True).stdout


def aligned_blocks(position: int, cigar: str) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    blocks = []
    introns = []
    cursor = position
    for length_text, op in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op in {"M", "=", "X"}:
            blocks.append((cursor, cursor + length - 1))
            cursor += length
        elif op == "N":
            introns.append((cursor, cursor + length - 1))
            cursor += length
        elif op == "D":
            cursor += length
    return blocks, introns


def overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    loci = parse_loci(args.loci)
    with args.samples.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if len(samples) != 9:
        raise ValueError(f"Expected nine RNA-seq samples, found {len(samples)}")

    per_sample: list[dict[str, object]] = []
    depth_vectors: dict[tuple[str, str], list[int]] = {}

    for sample in samples:
        run = sample["run"]
        bam = args.bam_root / run / f"{run}.backfill_targets.bam"
        if not bam.exists() or not Path(str(bam) + ".bai").exists():
            raise FileNotFoundError(f"Missing BAM or index: {bam}")
        for locus in loci:
            start = int(locus["coding_lo"])
            end = int(locus["coding_hi"])
            length = end - start + 1
            region = f"{locus['chrom']}:{start}-{end}"
            sam = run_text([args.samtools, "view", str(bam), region])
            fragments: set[str] = set()
            primary_alignments = 0
            splice_through = 0
            for line in sam.splitlines():
                fields = line.split("\t")
                if len(fields) < 11:
                    continue
                flag = int(fields[1])
                mapq = int(fields[4])
                if flag & 0x4 or flag & 0x100 or flag & 0x800 or mapq < 20:
                    continue
                blocks, introns = aligned_blocks(int(fields[3]), fields[5])
                if not any(overlap(a, b, start, end) >= 10 for a, b in blocks):
                    continue
                primary_alignments += 1
                fragments.add(fields[0])
                if any(overlap(a, b, start, end) > 0 for a, b in introns):
                    splice_through += 1

            depth_text = run_text([args.samtools, "depth", "-aa", "-r", region, str(bam)])
            depths = [int(line.split("\t")[2]) for line in depth_text.splitlines() if line]
            if len(depths) != length:
                raise ValueError(f"Depth vector length mismatch for {run} {region}: {len(depths)} vs {length}")
            depth_vectors[(run, str(locus["member_id"]))] = depths
            breadth1 = sum(value >= 1 for value in depths) / length
            breadth5 = sum(value >= 5 for value in depths) / length
            per_sample.append(
                {
                    "run": run,
                    "treatment": sample.get("treatment", sample.get("condition", "")),
                    "replicate": sample["replicate"],
                    "member_id": locus["member_id"],
                    "family": locus["family"],
                    "coding_interval": region,
                    "coding_length_bp": length,
                    "primary_mapq20_alignments": primary_alignments,
                    "unique_primary_mapq20_fragments": len(fragments),
                    "mean_depth": f"{sum(depths) / length:.6f}",
                    "maximum_depth": max(depths, default=0),
                    "breadth_depth_ge1": f"{breadth1:.6f}",
                    "breadth_depth_ge5": f"{breadth5:.6f}",
                    "primary_alignments_spliced_through_coding_interval": splice_through,
                }
            )

    locus_rows: list[dict[str, object]] = []
    for locus in loci:
        member_id = str(locus["member_id"])
        sample_rows = [row for row in per_sample if row["member_id"] == member_id]
        length = int(sample_rows[0]["coding_length_bp"])
        pooled_depth = [0] * length
        for sample in samples:
            vector = depth_vectors[(sample["run"], member_id)]
            pooled_depth = [a + b for a, b in zip(pooled_depth, vector)]
        supporting_samples = sum(
            int(row["unique_primary_mapq20_fragments"]) >= 3 and float(row["breadth_depth_ge1"]) >= 0.5
            for row in sample_rows
        )
        total_fragments = sum(int(row["unique_primary_mapq20_fragments"]) for row in sample_rows)
        pooled_breadth1 = sum(value >= 1 for value in pooled_depth) / length
        pooled_breadth5 = sum(value >= 5 for value in pooled_depth) / length
        if supporting_samples >= 2 and pooled_breadth1 >= 0.8:
            evidence = "CONSISTENT_SHORT_READ_LOCUS_SUPPORT"
        elif total_fragments >= 3 and pooled_breadth1 >= 0.5:
            evidence = "SPARSE_OR_SINGLE_SAMPLE_SHORT_READ_SUPPORT"
        elif total_fragments > 0:
            evidence = "WEAK_LOCAL_OVERLAP_ONLY"
        else:
            evidence = "NO_SHORT_READ_SUPPORT_DETECTED"

        condition_fragments = defaultdict(int)
        for row in sample_rows:
            condition_fragments[str(row["treatment"])] += int(row["unique_primary_mapq20_fragments"])
        locus_rows.append(
            {
                "member_id": member_id,
                "family": locus["family"],
                "historical_locus": locus["historical_locus"],
                "coding_interval": sample_rows[0]["coding_interval"],
                "total_unique_fragment_counts_summed_across_samples": total_fragments,
                "samples_with_fragment_and_ge50pct_breadth_support": supporting_samples,
                "pooled_mean_depth": f"{sum(pooled_depth) / length:.6f}",
                "pooled_breadth_depth_ge1": f"{pooled_breadth1:.6f}",
                "pooled_breadth_depth_ge5": f"{pooled_breadth5:.6f}",
                "CK_fragment_count_sum": condition_fragments["CK"],
                "S200_fragment_count_sum": condition_fragments["S200"],
                "S400_fragment_count_sum": condition_fragments["S400"],
                "short_read_evidence_class": evidence,
                "interpretation_limit": "Short-read overlap does not by itself establish a complete ORF, transcript model, or functional AMP gene",
            }
        )

    write_tsv(args.per_sample_output, per_sample)
    write_tsv(args.locus_output, locus_rows)
    print(f"Wrote {len(per_sample)} per-sample rows and {len(locus_rows)} locus summaries")


if __name__ == "__main__":
    main()
