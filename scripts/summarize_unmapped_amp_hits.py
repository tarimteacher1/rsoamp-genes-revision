#!/usr/bin/env python3
"""Summarize strict targeted AMP similarities without calling novel genes."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fields = [
        "run", "fragment_id", "mate", "subject_id", "subject_class", "identity_pct",
        "alignment_length_aa", "query_length_nt", "subject_length_aa", "query_coverage_pct",
        "subject_coverage_pct", "evalue", "bitscore", "interpretation",
    ]
    rows = []
    for sample in sorted(path for path in args.classification_root.iterdir() if path.is_dir()):
        path = sample / "residual_genome_unmapped_vs_AMP.diamond.tsv"
        if not path.exists():
            raise FileNotFoundError(path)
        best = {}
        with path.open(encoding="utf-8") as handle:
            for values in csv.reader(handle, delimiter="\t"):
                qseqid, subject, pident, length, qlen, slen, qcov, scov, evalue, bitscore = values
                if float(pident) < 50 or int(length) < 20 or float(qcov) < 50 or float(scov) < 30:
                    continue
                if qseqid not in best or float(bitscore) > float(best[qseqid][-1]):
                    best[qseqid] = values
        for qseqid, values in best.items():
            qseqid, subject, pident, length, qlen, slen, qcov, scov, evalue, bitscore = values
            fragment, mate = qseqid.rsplit("/", 1)
            subject_class = "audited_Rso_primary_catalogue" if subject.startswith("RSO_PRIMARY|") else "reviewed_cross_species_reference"
            interpretation = (
                "known_Rso_AMP_similarity_in_genome_unmapped_read; investigate reference divergence or read quality"
                if subject_class == "audited_Rso_primary_catalogue"
                else "curated_AMP_similarity_only; isolated read evidence is insufficient for a novel gene call"
            )
            rows.append(
                dict(
                    zip(
                        fields,
                        [sample.name, fragment, mate, subject, subject_class, pident, length, qlen, slen, qcov, scov, evalue, bitscore, interpretation],
                    )
                )
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["run"], row["fragment_id"], row["mate"])))
    known = sum(row["subject_class"] == "audited_Rso_primary_catalogue" for row in rows)
    cross = len(rows) - known
    args.summary.write_text(
        "metric\tvalue\n"
        f"strict_read_level_AMP_hits\t{len(rows)}\n"
        f"hits_to_audited_Rso_primary_catalogue\t{known}\n"
        f"hits_to_reviewed_cross_species_AMP_reference\t{cross}\n"
        "high_confidence_novel_AMP_gene_calls\t0\n"
        "novel_gene_call_rule\tRead-level similarity alone was not treated as a complete ORF or gene model\n",
        encoding="utf-8",
    )
    print(f"PASS strict AMP read hits={len(rows)} novel_gene_calls=0")


if __name__ == "__main__":
    main()
