#!/usr/bin/env python3
"""Prepare complete start-to-stop peptides from recursive genome-screen ORFs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ORF_COORD_RE = re.compile(r"\[(\d+)\s+-\s+(\d+)\]")


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    identifier = ""
    description = ""
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier:
                    records.append((identifier, description, "".join(parts).upper()))
                description = line[1:]
                identifier = description.split()[0]
                parts = []
            else:
                parts.append(line)
    if identifier:
        records.append((identifier, description, "".join(parts).upper()))
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_domains(path: Path, accession: str) -> set[str]:
    identifiers: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.split(maxsplit=22)
            if len(fields) >= 5 and fields[4] == accession:
                identifiers.add(fields[0])
    return identifiers


def read_windows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for index, raw in enumerate(handle, 1):
            chrom, start0, end0, name = raw.rstrip("\n").split("\t")[:4]
            cluster_id, expected_strand = name.split("|", 1)
            rows.append(
                {
                    "index": index,
                    "cluster_id": cluster_id,
                    "chrom": chrom,
                    "start0": int(start0),
                    "end0": int(end0),
                    "expected_strand": expected_strand,
                }
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows-bed", required=True, type=Path)
    parser.add_argument("--windows-fasta", required=True, type=Path)
    parser.add_argument("--orfs-fasta", required=True, type=Path)
    parser.add_argument("--domains", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--domain-accession", default="PF00234.28")
    args = parser.parse_args()

    windows = read_windows(args.windows_bed)
    window_sequences = read_fasta(args.windows_fasta)
    if len(windows) != len(window_sequences):
        raise ValueError("Window BED and FASTA counts differ")
    domain_orfs = read_domains(args.domains, args.domain_accession)
    orfs = {identifier: (description, sequence) for identifier, description, sequence in read_fasta(args.orfs_fasta)}
    if not domain_orfs <= set(orfs):
        raise ValueError("Domain-hit ORFs are missing from the ORF FASTA")

    rows: list[dict[str, object]] = []
    fasta_records: list[tuple[str, str]] = []
    for identifier in sorted(domain_orfs):
        prefix = identifier.split("_", 2)
        if len(prefix) < 3 or not prefix[1].isdigit():
            raise ValueError(f"Cannot identify source window for {identifier}")
        window_index = int(prefix[1])
        window = windows[window_index - 1]
        window_sequence = window_sequences[window_index - 1][2]
        description, peptide = orfs[identifier]
        match = ORF_COORD_RE.search(description)
        if not match:
            raise ValueError(f"Cannot parse ORF coordinates: {description}")
        raw_start, raw_end = map(int, match.groups())
        reverse = "REVERSE SENSE" in description
        strand = "-" if reverse else "+"
        if strand != window["expected_strand"]:
            strand_note = "OPPOSITE_TO_RECURSIVE_HSP"
        else:
            strand_note = "MATCHES_RECURSIVE_HSP"
        local_low, local_high = sorted((raw_start, raw_end))
        if local_high - local_low + 1 != len(peptide) * 3:
            raise ValueError(f"Coordinate/peptide length mismatch for {identifier}")
        first_m = peptide.find("M")
        if reverse:
            downstream_nt = reverse_complement(window_sequence[max(0, local_low - 4) : local_low - 1])
            local_cds_low = local_low
            local_cds_high = local_high - max(first_m, 0) * 3
        else:
            downstream_nt = window_sequence[local_high : local_high + 3]
            local_cds_low = local_low + max(first_m, 0) * 3
            local_cds_high = local_high
        downstream_stop = downstream_nt if downstream_nt in {"TAA", "TAG", "TGA"} else ""
        complete = first_m >= 0 and bool(downstream_stop)
        start_stop = peptide[first_m:] if complete else ""
        global_start = int(window["start0"]) + local_cds_low
        global_end = int(window["start0"]) + local_cds_high
        output_id = f'{window["cluster_id"]}_{identifier}|start_stop'
        if complete:
            fasta_records.append((output_id, start_stop))
        rows.append(
            {
                "recursive_cluster_id": window["cluster_id"],
                "source_orf_id": identifier,
                "output_id": output_id.removesuffix("|start_stop"),
                "chrom": window["chrom"],
                "genomic_start": global_start,
                "genomic_end": global_end,
                "strand": strand,
                "strand_consistency": strand_note,
                "stop_delimited_length_aa": len(peptide),
                "first_methionine_position_aa": first_m + 1 if first_m >= 0 else "NA",
                "downstream_stop_codon": downstream_stop or "ABSENT_OR_WINDOW_BOUNDARY",
                "complete_start_stop_suborf": "YES" if complete else "NO",
                "start_stop_length_aa": len(start_stop) if complete else "NA",
                "start_stop_cysteine_count": start_stop.count("C") if complete else "NA",
                "family_HMM_accession": args.domain_accession,
                "catalogue_action_before_family_review": "HOLD_OUTSIDE_PRIMARY_CATALOGUE",
            }
        )

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in fasta_records:
            handle.write(f">{identifier}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")
    with args.manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Prepared {len(fasta_records)} complete start-to-stop peptides from {len(rows)} domain-hit ORFs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
