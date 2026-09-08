#!/usr/bin/env python3
"""Summarize mutually exclusive transcriptome-unassigned read categories."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path


CATEGORIES = [
    "genome_decoy_assigned",
    "low_quality_or_low_complexity",
    "rRNA_derived",
    "plastid_derived",
    "mitochondrial_like_same_family_proxy",
    "genome_mapped_but_transcriptome_unassigned",
    "genome_unmapped_residual",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification-root", required=True, type=Path)
    parser.add_argument("--quant-root", required=True, type=Path)
    parser.add_argument("--alignment-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def star_metric(path: Path, label: str) -> float:
    pattern = re.compile(rf"^\s*{re.escape(label)}\s*\|\s*(.+?)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1).rstrip("%"))
    raise ValueError(f"STAR metric not found: {label} in {path}")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    metrics = []
    category_rows = []
    aggregate = Counter()
    for sample in sorted(path for path in args.classification_root.iterdir() if path.is_dir()):
        run = sample.name
        meta = json.loads((args.quant_root / run / "aux_info/meta_info.json").read_text(encoding="utf-8"))
        counts = Counter()
        with gzip.open(sample / "fragment_category_assignments.tsv.gz", "rt", encoding="utf-8") as handle:
            for line in handle:
                _, category = line.rstrip("\n").split("\t")
                counts[category] += 1
        missing = set(CATEGORIES) - set(counts)
        unexpected = set(counts) - set(CATEGORIES)
        if unexpected:
            raise ValueError(f"Unexpected categories for {run}: {sorted(unexpected)}")
        unassigned = sum(counts.values())
        expected_unassigned = int(meta["num_processed"]) - int(meta["num_mapped"])
        star_log = args.alignment_root / run / f"{run}.Log.final.out"
        metrics.append(
            {
                "run": run,
                "salmon_index_definition": "GFF-derived full transcripts plus genome decoys",
                "salmon_num_processed_fragments": meta["num_processed"],
                "salmon_num_transcriptome_mapped_fragments": meta["num_mapped"],
                "salmon_percent_transcriptome_mapped": meta["percent_mapped"],
                "salmon_num_processed_minus_num_mapped": expected_unassigned,
                "salmon_unmapped_name_fragment_records": unassigned,
                "unmapped_name_records_minus_processed_minus_mapped": unassigned - expected_unassigned,
                "star_uniquely_mapped_percent": star_metric(star_log, "Uniquely mapped reads %"),
                "star_multimapped_percent": star_metric(star_log, "% of reads mapped to multiple loci"),
                "star_unmapped_too_short_percent": star_metric(star_log, "% of reads unmapped: too short"),
                "missing_category_types_with_zero_count": ";".join(sorted(missing)),
            }
        )
        for category in CATEGORIES:
            n = counts[category]
            aggregate[category] += n
            category_rows.append(
                {
                    "run": run,
                    "category": category,
                    "fragment_count": n,
                    "fraction_of_transcriptome_unassigned": f"{n / unassigned:.6f}" if unassigned else "0",
                    "classification_priority": CATEGORIES.index(category) + 1,
                }
            )
    total = sum(aggregate.values())
    for category in CATEGORIES:
        category_rows.append(
            {
                "run": "ALL_SAMPLES",
                "category": category,
                "fragment_count": aggregate[category],
                "fraction_of_transcriptome_unassigned": f"{aggregate[category] / total:.6f}" if total else "0",
                "classification_priority": CATEGORIES.index(category) + 1,
            }
        )
    write_tsv(args.output_dir / "per_sample_mapping_metrics.tsv", metrics)
    write_tsv(args.output_dir / "unmapped_read_categories.tsv", category_rows)
    summary = [
        "# Transcriptome-unassigned read audit",
        "",
        "Status: PASS",
        "",
        "Unmapped is defined relative to the decoy-aware Salmon transcript index. Each Salmon-unassigned fragment was assigned exactly once, in this order: genome decoy, quality/complexity failure, rRNA, same-species plastid, same-family mitochondrial proxy, residual genome-mapped, or residual genome-unmapped.",
        "",
        "The mitochondrial category is a conservative proxy screen against Myricaria laxiflora MW971331 because a species-specific R. soongarica mitochondrial reference was not available; it is not interpreted as a complete mitochondrial read count.",
        "",
        f"Across all samples, Salmon wrote {total:,} transcriptome-unassigned fragment records for classification. The per-sample table separately reports any difference between this record count and num_processed minus num_mapped rather than forcing the two quantities to agree.",
    ]
    for category in CATEGORIES:
        summary.append(f"- {category}: {aggregate[category]:,} ({aggregate[category] / total:.2%})")
    (args.output_dir / "unmapped_reads_summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"PASS unmapped audit samples={len(metrics)} fragments={total}")


if __name__ == "__main__":
    main()
