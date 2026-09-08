#!/usr/bin/env python3
"""Summarize primary long-read alignment support for six historical AMP ORFs."""

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
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--samtools", required=True)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def reference_blocks(pos_1based: int, cigar: str) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Return aligned query/reference blocks and skipped-reference introns, 1-based inclusive."""
    cursor = pos_1based
    blocks: list[tuple[int, int]] = []
    introns: list[tuple[int, int]] = []
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
        elif op in {"I", "S", "H", "P"}:
            continue
    return blocks, introns


def overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def load_loci(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 6:
        raise ValueError(f"Expected six loci, found {len(rows)}")
    for row in rows:
        chrom, remainder = row["historical_locus"].split(":", 1)
        match = re.match(r"(\d+)-(\d+)\(([+-])\)", remainder)
        if not match:
            raise ValueError(f"Cannot parse locus: {row['historical_locus']}")
        row["chrom"] = chrom
        row["strand"] = match.group(3)
    return rows


def sam_rows(samtools: str, bam: Path, region: str) -> list[list[str]]:
    result = subprocess.run(
        [samtools, "view", str(bam), region],
        check=True,
        text=True,
        capture_output=True,
    )
    return [line.split("\t") for line in result.stdout.splitlines() if line]


def main() -> None:
    args = parse_args()
    loci = load_loci(args.loci)
    output_rows: list[dict[str, object]] = []

    for locus in loci:
        chrom = locus["chrom"]
        coding_start = min(int(locus["coding_start"]), int(locus["coding_end"]))
        coding_end = max(int(locus["coding_start"]), int(locus["coding_end"]))
        coding_len = coding_end - coding_start + 1
        window_start = max(1, coding_start - 10000)
        window_end = coding_end + 10000
        region = f"{chrom}:{window_start}-{window_end}"

        metrics = defaultdict(int)
        max_fraction = 0.0
        supporting_reads: list[str] = []
        seen_primary: set[str] = set()

        for fields in sam_rows(args.samtools, args.bam, region):
            if len(fields) < 11:
                continue
            read_id = fields[0]
            flag = int(fields[1])
            mapq = int(fields[4])
            cigar = fields[5]
            if flag & 0x4 or flag & 0x100 or flag & 0x800 or mapq < 20:
                continue
            if read_id in seen_primary:
                continue
            seen_primary.add(read_id)
            metrics["primary_mapq20_window"] += 1

            read_strand = "-" if flag & 0x10 else "+"
            strand_match = read_strand == locus["strand"]
            if strand_match:
                metrics["strand_match"] += 1

            blocks, introns = reference_blocks(int(fields[3]), cigar)
            aligned_overlap = sum(overlap(a, b, coding_start, coding_end) for a, b in blocks)
            aligned_overlap = min(aligned_overlap, coding_len)
            fraction = aligned_overlap / coding_len
            max_fraction = max(max_fraction, fraction)
            if aligned_overlap > 0:
                metrics["overlap_coding"] += 1
            if strand_match and fraction >= 0.5:
                metrics["strand_cover50"] += 1
            if strand_match and fraction >= 0.8:
                metrics["strand_cover80"] += 1
                if len(supporting_reads) < 20:
                    supporting_reads.append(read_id)

            max_contiguous = max(
                (overlap(a, b, coding_start, coding_end) for a, b in blocks),
                default=0,
            )
            if strand_match and max_contiguous / coding_len >= 0.8:
                metrics["strand_contiguous80"] += 1
            if strand_match and any(
                intron_start <= coding_end and intron_end >= coding_start
                for intron_start, intron_end in introns
            ):
                metrics["splice_through_coding"] += 1

        if metrics["strand_contiguous80"] >= 2:
            evidence = "STRONG_CONTIGUOUS_TRANSCRIPT_SUPPORT"
        elif metrics["strand_contiguous80"] == 1 or metrics["strand_cover50"] >= 2:
            evidence = "MODERATE_TRANSCRIPT_SUPPORT"
        elif metrics["overlap_coding"] >= 1:
            evidence = "WEAK_OR_STRAND_INCONSISTENT_SUPPORT"
        else:
            evidence = "NO_TARGET_OVERLAP_DETECTED"

        output_rows.append(
            {
                "member_id": locus["member_id"],
                "family": locus["family"],
                "locus": locus["historical_locus"],
                "coding_interval": f"{chrom}:{coding_start}-{coding_end}({locus['strand']})",
                "coding_length_bp": coding_len,
                "primary_mapq20_reads_in_10kb_window": metrics["primary_mapq20_window"],
                "strand_matching_reads_in_10kb_window": metrics["strand_match"],
                "primary_reads_overlapping_coding_interval": metrics["overlap_coding"],
                "strand_matching_reads_covering_ge50pct": metrics["strand_cover50"],
                "strand_matching_reads_covering_ge80pct": metrics["strand_cover80"],
                "strand_matching_contiguous_reads_covering_ge80pct": metrics["strand_contiguous80"],
                "strand_matching_reads_spliced_through_coding_interval": metrics["splice_through_coding"],
                "maximum_aligned_coding_fraction": f"{max_fraction:.4f}",
                "supporting_read_ids": ";".join(supporting_reads),
                "transcript_evidence_class": evidence,
                "interpretation_limit": "Alignment support alone does not establish a complete ORF or functional AMP gene",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()
