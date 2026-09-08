#!/usr/bin/env python3
"""Build evidence-bounded qPCR and reference-gene candidate packets."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


RUN_GROUPS = {
    "CK": ["SRR27540881", "SRR27540880", "SRR27540879"],
    "S200": ["SRR27540878", "SRR27540877", "SRR27540876"],
    "S400": ["SRR27540875", "SRR27540884", "SRR27540883"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--de-results", required=True, type=Path)
    parser.add_argument("--amp-tpm", required=True, type=Path)
    parser.add_argument("--mappability", required=True, type=Path)
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--all-normalized-counts", required=True, type=Path)
    parser.add_argument("--all-de-results", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line)
    return {key: "".join(chunks) for key, chunks in records.items()}


def finite_float(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def true_value(value: str) -> bool:
    return value.strip().upper() in {"TRUE", "T", "YES", "1"}


def means_by_group(row: dict[str, str]) -> dict[str, float | None]:
    result: dict[str, float | None] = {}
    for group, runs in RUN_GROUPS.items():
        values = [finite_float(row.get(run, "")) for run in runs]
        observed = [value for value in values if value is not None]
        result[group] = statistics.mean(observed) if observed else None
    return result


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalogue = read_table(args.catalogue, "\t")
    de_rows = read_table(args.de_results, "\t")
    tpm_rows = read_table(args.amp_tpm, "\t")
    mappability_rows = read_table(args.mappability, "\t")
    transcripts = read_fasta(args.transcripts)

    catalogue_by_member = {row["member_id"]: row for row in catalogue}
    de_by_member_contrast = {(row["member_id"], row["contrast"]): row for row in de_rows}
    tpm_by_member = {row["member_id"]: row for row in tpm_rows}
    mappability_by_member = {row["member_id"]: row for row in mappability_rows}

    candidates = []
    for member_id, member in catalogue_by_member.items():
        s400 = de_by_member_contrast.get((member_id, "S400_vs_CK"))
        if not s400 or not true_value(s400.get("DE_call_FDR_lt_0.05_abs_log2FC_gt_1", "")):
            continue
        fold_change = finite_float(s400.get("log2FoldChange", ""))
        adjusted_p = finite_float(s400.get("padj_BH", ""))
        if fold_change is None or adjusted_p is None:
            continue
        mapping = mappability_by_member.get(member_id, {})
        risk = mapping.get("mappability_risk", "NOT_ASSESSED")
        risk_rank = {"LOW": 0, "MODERATE": 1, "HIGH": 2}.get(risk, 3)
        candidates.append(
            {
                "member_id": member_id,
                "family": member["family"],
                "direction": "UP" if fold_change > 0 else "DOWN",
                "risk_rank": risk_rank,
                "padj": adjusted_p,
                "abs_log2fc": abs(fold_change),
            }
        )
    candidates.sort(key=lambda row: (row["risk_rank"], row["padj"], -row["abs_log2fc"]))

    quotas = [
        ("nsLTP", "UP", 3),
        ("Defensin", "UP", 2),
        ("Snakin_GASA", "DOWN", 2),
    ]
    selected: list[dict[str, object]] = []
    selected_ids: set[str] = set()
    for family, direction, quota in quotas:
        eligible = [
            row for row in candidates
            if row["family"] == family and row["direction"] == direction and row["member_id"] not in selected_ids
        ]
        for row in eligible[:quota]:
            selected.append(row)
            selected_ids.add(str(row["member_id"]))
    for row in candidates:
        if len(selected) >= 8:
            break
        if str(row["member_id"]) not in selected_ids:
            selected.append(row)
            selected_ids.add(str(row["member_id"]))
    if len(selected) < min(6, len(candidates)):
        raise ValueError("Could not construct the required family-balanced qPCR candidate panel")

    packet_rows: list[dict[str, object]] = []
    sequence_rows = []
    for rank, selection in enumerate(selected, start=1):
        member_id = str(selection["member_id"])
        member = catalogue_by_member[member_id]
        s200 = de_by_member_contrast.get((member_id, "S200_vs_CK"), {})
        s400 = de_by_member_contrast[(member_id, "S400_vs_CK")]
        mapping = mappability_by_member.get(member_id, {})
        tpm = tpm_by_member.get(member_id, {})
        group_means = means_by_group(tpm)
        near_identical = int(mapping.get("near_identical_nonself_hits_95pct_80pct_coverage", "0") or 0)
        risk = mapping.get("mappability_risk", "NOT_ASSESSED")
        unique_region = (
            "LIKELY_AVAILABLE_TRANSCRIPT_LEVEL; PRIMER_BLAST_REQUIRED"
            if risk == "LOW" and near_identical == 0
            else "UNCERTAIN; TARGET_SPECIFIC_PRIMER_DESIGN_AND_IN_SILICO_PCR_REQUIRED"
        )
        rationale = (
            f"{selection['direction'].lower()}regulated S400-versus-CK transcriptome candidate; "
            f"family-balanced selection; transcript-level mappability risk={risk.lower()}"
        )
        packet_rows.append(
            {
                "proposed_priority": rank,
                "member_id": member_id,
                "family": member["family"],
                "canonical_status": member["final_subclass"],
                "gene_model_confidence": member["classification_confidence"],
                "gene_id": member["gene_id"],
                "transcript_id": member["protein_id"],
                "S200_log2FC": s200.get("log2FoldChange", ""),
                "S200_padj_BH": s200.get("padj_BH", ""),
                "S400_log2FC": s400.get("log2FoldChange", ""),
                "S400_padj_BH": s400.get("padj_BH", ""),
                "mean_TPM_CK": f"{group_means['CK']:.6f}" if group_means["CK"] is not None else "",
                "mean_TPM_S200": f"{group_means['S200']:.6f}" if group_means["S200"] is not None else "",
                "mean_TPM_S400": f"{group_means['S400']:.6f}" if group_means["S400"] is not None else "",
                "mappability_risk": risk,
                "best_nonself_identity_pct": mapping.get("best_nonself_identity_pct", ""),
                "best_nonself_query_coverage": mapping.get("best_nonself_query_coverage", ""),
                "median_salmon_ambiguous_fraction": mapping.get("median_salmon_ambiguous_fraction", ""),
                "unique_region_availability": unique_region,
                "selection_rationale": rationale,
                "validation_status": "USER-WETLAB-PENDING; NOT_EXPERIMENTALLY_VALIDATED",
            }
        )
        transcript_id = member["protein_id"]
        sequence = transcripts.get(transcript_id)
        if not sequence:
            raise ValueError(f"Selected qPCR candidate lacks transcript sequence: {member_id} {transcript_id}")
        sequence_rows.append((member_id, transcript_id, sequence))

    write_tsv(args.output_dir / "qPCR_candidate_packet.tsv", packet_rows)
    with (args.output_dir / "qPCR_candidate_sequences.fasta").open("w", encoding="ascii") as handle:
        for member_id, transcript_id, sequence in sequence_rows:
            handle.write(f">{member_id} transcript_id={transcript_id}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    normalized = read_table(args.all_normalized_counts, ",")
    all_de = read_table(args.all_de_results, ",")
    de_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in all_de:
        if row.get("contrast"):
            de_by_gene[row["gene_id"]].append(row)
    reference_rows = []
    for row in normalized:
        gene_id = row.get("gene_id") or next(iter(row.values()))
        values = [finite_float(row.get(run, "")) for runs in RUN_GROUPS.values() for run in runs]
        observed = [value for value in values if value is not None]
        if len(observed) != 9 or statistics.mean(observed) < 100:
            continue
        mean_value = statistics.mean(observed)
        cv = statistics.stdev(observed) / mean_value if mean_value else float("inf")
        group_means = means_by_group(row)
        group_values = [value for value in group_means.values() if value is not None]
        if len(group_values) != 3 or min(group_values) <= 0:
            continue
        max_condition_ratio = max(group_values) / min(group_values)
        tests = de_by_gene.get(gene_id, [])
        stable_tests = bool(tests) and all(
            (finite_float(test.get("padj", "")) is None or float(test["padj"]) >= 0.5)
            and (finite_float(test.get("log2FoldChange", "")) is None or abs(float(test["log2FoldChange"])) <= 0.25)
            for test in tests
        )
        if cv <= 0.15 and max_condition_ratio <= 1.20 and stable_tests:
            reference_rows.append(
                {
                    "gene_id": gene_id,
                    "mean_normalized_count": f"{mean_value:.6f}",
                    "sample_CV": f"{cv:.6f}",
                    "max_to_min_condition_mean_ratio": f"{max_condition_ratio:.6f}",
                    "all_three_DE_contrasts_padj_ge_0.5_or_NA": "YES",
                    "all_three_abs_log2FC_le_0.25_or_NA": "YES",
                    "candidate_label": "RNA-seq-based candidate reference gene requiring experimental stability validation",
                    "validation_status": "USER-WETLAB-PENDING",
                }
            )
    reference_rows.sort(key=lambda row: (float(row["sample_CV"]), -float(row["mean_normalized_count"])))
    if not reference_rows:
        raise ValueError("No reference-gene candidate passed the predeclared RNA-seq stability screen")
    write_tsv(args.output_dir / "reference_gene_candidates.tsv", reference_rows[:12])

    (args.output_dir / "QPCR_CANDIDATE_PACKET_COMPLETE.PASS").write_text(
        f"PASS\tqPCR_candidates={len(packet_rows)}\treference_gene_candidates={min(12, len(reference_rows))}\n",
        encoding="utf-8",
    )
    print(f"PASS qPCR candidates={len(packet_rows)} reference candidates={min(12, len(reference_rows))}")


if __name__ == "__main__":
    main()
