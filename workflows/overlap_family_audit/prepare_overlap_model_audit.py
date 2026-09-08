#!/usr/bin/env python3
"""Prepare non-intergenic M6 transcript models for family-boundary review.

The focused discovery output distinguishes intergenic models from models that
overlap or alter an existing reference transcript.  The latter cannot create a
new locus count without additional review, but complete start-to-stop ORFs may
still expose an alternative isoform or a family-boundary problem.  This script
extracts those ORFs without changing the frozen catalogue.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                records[current].append(line)
    return {identifier: "".join(parts) for identifier, parts in records.items()}


def write_fasta(path: Path, rows: list[tuple[str, str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, description, sequence in rows:
            handle.write(f">{identifier} {description}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-results", required=True, type=Path)
    parser.add_argument("--catalogue-all", required=True, type=Path)
    parser.add_argument("--nsltp-evidence", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    source = args.source_results.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    evidence = [
        row
        for row in read_tsv(source / "candidate_evidence.tsv")
        if row["model_origin"] == "OTHER_OR_ALTERED_REFERENCE_RELATIONSHIP"
    ]
    if len(evidence) != 32:
        raise ValueError(f"Expected 32 overlap/altered evidence rows, observed {len(evidence)}")
    proteins = read_fasta(source / "screen/stop_delimited_orfs.faa")
    missing = sorted({row["orf_id"] for row in evidence} - proteins.keys())
    if missing:
        raise ValueError(f"Candidate ORFs missing from source FASTA: {missing[:5]}")

    audited_by_protein = {
        row["protein_id"]: row for row in read_tsv(args.catalogue_all.resolve()) if row.get("protein_id")
    }
    nsltp_by_protein = {
        row["protein_id"]: row for row in read_tsv(args.nsltp_evidence.resolve()) if row.get("protein_id")
    }
    all_orfs: list[tuple[str, str, str]] = []
    start_stop: list[tuple[str, str, str]] = []
    manifest: list[dict[str, object]] = []
    for index, row in enumerate(
        sorted(evidence, key=lambda value: (value["chrom"], int(value["start"]), value["transcript_id"], value["orf_id"])),
        start=1,
    ):
        oid = row["orf_id"]
        sequence = proteins[oid]
        first_m = int(row["first_methionine_position_aa"]) if row["first_methionine_position_aa"] != "NA" else 0
        complete = row["possible_start_stop_suborf"] == "True" and first_m > 0 and row["downstream_stop_present"] == "True"
        peptide = sequence[first_m - 1 :] if complete else ""
        existing = audited_by_protein.get(row["reference_id"], {})
        nsltp_existing = nsltp_by_protein.get(row["reference_id"], {})
        audit_id = f"R1M6_OVERLAP{index:02d}"
        all_orfs.append((oid, f"audit_id={audit_id};transcript_id={row['transcript_id']}", sequence))
        if complete:
            start_stop.append(
                (
                    oid + "|start_stop",
                    f"audit_id={audit_id};transcript_id={row['transcript_id']};reference_id={row['reference_id']}",
                    peptide,
                )
            )
        manifest.append(
            {
                "audit_id": audit_id,
                "transcript_id": row["transcript_id"],
                "orf_id": oid,
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "gffcompare_class": row["gffcompare_class"],
                "reference_gene_id": row["reference_gene_id"],
                "reference_id": row["reference_id"],
                "reference_in_audited_amp_pool": "YES" if existing else "NO",
                "existing_member_id": existing.get("member_id", ""),
                "existing_family": existing.get("family", ""),
                "existing_catalogue_status": existing.get("catalogue_status", ""),
                "reference_in_frozen_nsltp_pool": "YES" if nsltp_existing else "NO",
                "frozen_nsltp_final_classification": nsltp_existing.get("final_classification", ""),
                "frozen_nsltp_rule_trace": nsltp_existing.get("classification_rule_trace", ""),
                "review_priority": row["review_priority"],
                "focused_subset_reliable_samples": row["samples_with_3_unique_fragments_and_70pct_exon_coverage"],
                "stop_delimited_orf_length_aa": len(sequence),
                "complete_start_stop_suborf": "YES" if complete else "NO",
                "start_stop_length_aa": len(peptide) if peptide else "NA",
                "start_stop_cysteine_count": peptide.count("C") if peptide else "NA",
                "discovery_family_HMMs": row["family_HMMs"],
                "discovery_best_AMP_reference": row["best_AMP_reference"],
                "discovery_best_AMP_bitscore": row["best_AMP_bitscore"],
                "catalogue_action": "NO_AUTOMATIC_CHANGE",
            }
        )

    write_fasta(out / "overlap_all_stop_delimited_orfs.faa", all_orfs)
    write_fasta(out / "overlap_complete_start_stop_suborfs.faa", start_stop)
    fields = list(manifest[0])
    with (out / "overlap_model_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)

    summary = {
        "overlap_or_altered_ORF_rows": len(evidence),
        "complete_start_stop_ORFs": len(start_stop),
        "reference_relationships": dict(sorted(Counter(row["gffcompare_class"] for row in evidence).items())),
        "existing_audited_AMP_reference_rows": sum(row["reference_in_audited_amp_pool"] == "YES" for row in manifest),
        "existing_frozen_nsltp_pool_reference_rows": sum(row["reference_in_frozen_nsltp_pool"] == "YES" for row in manifest),
        "catalogue_change": "NONE_AUTOMATIC",
        "scope": "Alternative/overlapping transcript model and family-boundary audit; not a new-locus call or antimicrobial activity test",
    }
    (out / "overlap_model_preparation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
