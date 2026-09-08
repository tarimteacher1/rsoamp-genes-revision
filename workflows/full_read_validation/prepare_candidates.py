#!/usr/bin/env python3
"""Prepare an auditable intergenic-candidate set for R1.M6 validation."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


TRANSCRIPT_RE = re.compile(r'transcript_id "([^"]+)"')


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA ID: {current}")
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence before first FASTA header: {path}")
            else:
                records[current].append(line)
    return {key: "".join(value) for key, value in records.items()}


def write_fasta(path: Path, records: list[tuple[str, str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, description, sequence in records:
            handle.write(f">{identifier} {description}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def load_candidates(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row["model_origin"] == "INTERGENIC_CANDIDATE"]
    if not selected:
        raise ValueError("No intergenic candidates found")
    return selected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-results", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--window-bp", type=int, default=10_000)
    args = parser.parse_args()

    source = args.source_results.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    evidence_path = source / "candidate_evidence.tsv"
    gtf_path = source / "screen/comparison.annotated.gtf"
    orf_path = source / "screen/stop_delimited_orfs.faa"
    candidates = load_candidates(evidence_path)
    transcript_ids = {row["transcript_id"] for row in candidates}
    priority_ids = {
        row["orf_id"]
        for row in candidates
        if row["review_priority"] == "FULL_READ_MODEL_VALIDATION_PRIORITY"
    }

    selected_gtf: list[str] = []
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    transcript_rows: set[str] = set()
    with gtf_path.open(encoding="utf-8") as handle:
        for raw in handle:
            match = TRANSCRIPT_RE.search(raw)
            if not match or match.group(1) not in transcript_ids:
                continue
            selected_gtf.append(raw)
            fields = raw.rstrip("\n").split("\t")
            tid = match.group(1)
            if fields[2] == "transcript":
                transcript_rows.add(tid)
            elif fields[2] == "exon":
                exons[tid].append((int(fields[3]), int(fields[4])))
    if transcript_rows != transcript_ids or set(exons) != transcript_ids:
        raise ValueError("Candidate transcript/exon identities do not join to the GTF")
    (out / "intergenic_candidate_models.gtf").write_text("".join(selected_gtf), encoding="utf-8")

    proteins = read_fasta(orf_path)
    if not {row["orf_id"] for row in candidates} <= proteins.keys():
        raise ValueError("Candidate ORF identities do not join to the ORF FASTA")

    all_records: list[tuple[str, str, str]] = []
    start_stop_records: list[tuple[str, str, str]] = []
    manifest: list[dict[str, object]] = []
    for row in candidates:
        oid = row["orf_id"]
        sequence = proteins[oid]
        first_m = int(row["first_methionine_position_aa"]) if row["first_methionine_position_aa"] != "NA" else 0
        trimmed = sequence[first_m - 1 :] if first_m else ""
        all_records.append((oid, f'transcript_id={row["transcript_id"]};priority={row["review_priority"]}', sequence))
        if oid in priority_ids:
            if not first_m or row["downstream_stop_present"] != "True":
                raise ValueError(f"Priority ORF lacks the recorded start/stop indicators: {oid}")
            start_stop_records.append((oid + "|start_stop", f'transcript_id={row["transcript_id"]};source_orf={oid}', trimmed))
        cys_positions = [index for index, aa in enumerate(trimmed, start=1) if aa == "C"]
        first25 = trimmed[:25]
        hydrophobic = sum(aa in "AILMFWVY" for aa in first25) / len(first25) if first25 else 0.0
        manifest.append(
            {
                "orf_id": oid,
                "transcript_id": row["transcript_id"],
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "review_priority": row["review_priority"],
                "stop_delimited_length_aa": len(sequence),
                "first_methionine_position_aa": first_m or "NA",
                "start_stop_length_aa": len(trimmed) if trimmed else "NA",
                "start_stop_cysteine_count": len(cys_positions) if trimmed else "NA",
                "start_stop_cysteine_positions": ",".join(map(str, cys_positions)),
                "start_stop_first25_hydrophobic_fraction": f"{hydrophobic:.3f}" if first25 else "NA",
                "family_HMMs": row["family_HMMs"],
                "best_AMP_reference": row["best_AMP_reference"],
                "best_AMP_bitscore": row["best_AMP_bitscore"],
                "best_boundary_control": row["best_boundary_control"],
                "best_control_bitscore": row["best_control_bitscore"],
            }
        )

    write_fasta(out / "intergenic_stop_delimited_orfs.faa", all_records)
    write_fasta(out / "priority_start_stop_suborfs.faa", start_stop_records)
    fields = list(manifest[0])
    with (out / "candidate_sequence_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)

    with (out / "candidate_exact_exons.bed").open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(candidates, key=lambda value: (value["chrom"], int(value["start"]), value["transcript_id"])):
            for exon_number, (start, end) in enumerate(sorted(exons[row["transcript_id"]]), start=1):
                handle.write(
                    f'{row["chrom"]}\t{start - 1}\t{end}\t{row["transcript_id"]}.exon{exon_number}\t0\t{row["strand"]}\n'
                )

    intervals: list[tuple[str, int, int, set[str]]] = []
    for row in candidates:
        intervals.append(
            (
                row["chrom"],
                max(0, int(row["start"]) - 1 - args.window_bp),
                int(row["end"]) + args.window_bp,
                {row["transcript_id"]},
            )
        )
    merged: list[tuple[str, int, int, set[str]]] = []
    for chrom, start, end, ids in sorted(intervals, key=lambda value: (value[0], value[1], value[2])):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            old_chrom, old_start, old_end, old_ids = merged[-1]
            merged[-1] = (old_chrom, old_start, max(old_end, end), old_ids | ids)
        else:
            merged.append((chrom, start, end, set(ids)))
    with (out / "candidate_windows_10kb.bed").open("w", encoding="utf-8", newline="\n") as handle:
        for index, (chrom, start, end, ids) in enumerate(merged, start=1):
            handle.write(f'{chrom}\t{start}\t{end}\tM6WIN{index:02d}|{";".join(sorted(ids))}\n')

    summary = {
        "source_results": str(source),
        "intergenic_orf_rows": len(candidates),
        "intergenic_transcript_models": len(transcript_ids),
        "priority_start_stop_orfs": len(priority_ids),
        "merged_mapping_windows": len(merged),
        "window_extension_bp": args.window_bp,
        "catalogue_change": "NONE",
    }
    (out / "candidate_preparation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
