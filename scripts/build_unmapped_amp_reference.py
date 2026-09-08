#!/usr/bin/env python3
"""Build a labeled protein reference for targeted AMP searches in residual reads."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    return parser.parse_args()


def fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].split()[0]
            records[current] = []
        elif current is not None:
            records[current].append(line)
    return {key: "".join(value).replace("*", "") for key, value in records.items()}


def main() -> None:
    args = parse_args()
    records: list[dict[str, str]] = []
    primary = fasta(args.revision / "03_annotation_rescue/final_amp_primary_catalogue.faa")
    for member, sequence in primary.items():
        records.append(
            {
                "reference_id": f"RSO_PRIMARY|{member}",
                "reference_class": "audited_Rso_primary_catalogue",
                "source_id": member,
                "sequence": sequence,
            }
        )

    for family in ("Defensin", "Snakin_GASA"):
        path = args.revision / f"02_phylogeny/reference/{family}.uniprot_reviewed.filtered.tsv"
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                records.append(
                    {
                        "reference_id": f"UNIPROT_{family}|{row['Entry']}",
                        "reference_class": f"reviewed_{family}_reference",
                        "source_id": row["Entry"],
                        "sequence": row["Sequence"].replace("*", ""),
                    }
                )
    for entry, sequence in fasta(
        args.revision / "01_nslTP_curation/reference/positive_reviewed_arabidopsis_nsLTP.faa"
    ).items():
        records.append(
            {
                "reference_id": f"UNIPROT_nsLTP|{entry}",
                "reference_class": "reviewed_nsLTP_reference",
                "source_id": entry,
                "sequence": sequence,
            }
        )

    seen_sequences = set()
    unique = []
    for row in records:
        digest = hashlib.sha256(row["sequence"].encode()).hexdigest()
        if digest in seen_sequences:
            continue
        seen_sequences.add(digest)
        unique.append({**row, "sequence_sha256": digest})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii") as handle:
        for row in unique:
            handle.write(f">{row['reference_id']}\n{row['sequence']}\n")
    with args.manifest.open("w", encoding="utf-8", newline="") as handle:
        fields = ["reference_id", "reference_class", "source_id", "sequence_sha256"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({key: row[key] for key in fields} for row in unique)
    print(f"PASS targeted AMP references={len(unique)}")


if __name__ == "__main__":
    main()
