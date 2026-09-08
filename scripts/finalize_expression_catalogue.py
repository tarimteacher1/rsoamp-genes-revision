#!/usr/bin/env python3
"""Freeze expression outputs against the evidence-defined primary AMP catalogue."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--de-results", required=True, type=Path)
    parser.add_argument("--tpm", required=True, type=Path)
    parser.add_argument("--normalized-counts", required=True, type=Path)
    parser.add_argument("--vst", required=True, type=Path)
    parser.add_argument("--prefilter-audit", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
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


def de_call(value: str) -> bool:
    return value.strip().upper() in {"TRUE", "T", "1", "YES"}


def finite_number(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def test_status(audit: dict[str, str], result: dict[str, str] | None) -> str:
    retained = de_call(audit["retained_for_DE"])
    if retained != (result is not None):
        raise ValueError(f"DE row and prefilter audit disagree for {audit['gene_id']}")
    if not retained:
        reasons = {
            "FILTERED_ALL_ZERO_COUNTS": "NOT_TESTED_ALL_ZERO_COUNTS",
            "FILTERED_TOTAL_COUNT_LT_10": "NOT_TESTED_LOW_TOTAL_COUNT",
        }
        return reasons[audit["prefilter_status"]]
    if not finite_number(result["pvalue"]):
        return "P_VALUE_UNAVAILABLE"
    if not finite_number(result["padj"]):
        return "INDEPENDENTLY_FILTERED_PADJ_UNAVAILABLE"
    return "TESTED_WITH_ADJUSTED_P"


def matrix_subset(
    path: Path,
    catalogue: list[dict[str, str]],
    output: Path,
    audit: dict[str, dict[str, str]],
    prefiltered: bool,
) -> None:
    matrix = {row["gene_id"]: row for row in read_rows(path, ",")}
    value_columns = [column for column in next(iter(matrix.values())) if column != "gene_id"]
    rows = []
    for member in catalogue:
        values = matrix.get(member["gene_id"], {})
        status = audit[member["gene_id"]]
        expected = de_call(status["retained_for_DE"]) if prefiltered else True
        if bool(values) != expected:
            raise ValueError(f"Unexpected matrix membership in {path}: {member['gene_id']}")
        rows.append(
            {
                "member_id": member["member_id"],
                "family": member["family"],
                "final_subclass": member["final_subclass"],
                "gene_id": member["gene_id"],
                "quantification_status": "QUANTIFIED",
                "prefilter_status": status["prefilter_status"],
                "matrix_status": "AVAILABLE" if values else "NOT_AVAILABLE_PREFILTERED",
                **{column: values.get(column, "") for column in value_columns},
            }
        )
    write_tsv(output, rows)


def main() -> None:
    args = parse_args()
    catalogue = read_rows(args.catalogue, "\t")
    de_rows = read_rows(args.de_results, ",")
    if not catalogue:
        raise ValueError("Primary AMP catalogue is empty")
    de_lookup = {(row["gene_id"], row["contrast"]): row for row in de_rows}
    audit_rows = read_rows(args.prefilter_audit, "\t")
    audit = {row["gene_id"]: row for row in audit_rows}
    if len(audit) != len(audit_rows) or len(de_lookup) != len(de_rows):
        raise ValueError("Duplicate gene audit or gene/contrast result keys")
    if any(member["gene_id"] not in audit for member in catalogue):
        raise ValueError("A primary gene is absent from the prefilter audit")
    contrasts = ["S200_vs_CK", "S400_vs_CK", "S400_vs_S200"]

    final_rows: list[dict[str, object]] = []
    for member in catalogue:
        for contrast in contrasts:
            result = de_lookup.get((member["gene_id"], contrast))
            status = audit[member["gene_id"]]
            de_status = test_status(status, result)
            final_rows.append(
                {
                    "member_id": member["member_id"],
                    "family": member["family"],
                    "final_subclass": member["final_subclass"],
                    "final_subtype": member["final_subtype"],
                    "gene_id": member["gene_id"],
                    "protein_id": member["protein_id"],
                    "contrast": contrast,
                    "quantification_status": "QUANTIFIED",
                    "prefilter_status": status["prefilter_status"],
                    "DE_test_status": de_status,
                    "baseMean": result["baseMean"] if result else "",
                    "log2FoldChange": result["log2FoldChange"] if result else "",
                    "log2FC_apeglm": result["log2FC_apeglm"] if result else "",
                    "lfcSE": result["lfcSE"] if result else "",
                    "stat": result["stat"] if result else "",
                    "pvalue": result["pvalue"] if result else "",
                    "padj_BH": result["padj"] if result else "",
                    "DE_call_FDR_lt_0.05_abs_log2FC_gt_1": result["DE_call"] if de_status == "TESTED_WITH_ADJUSTED_P" else "NA",
                }
            )
    write_tsv(args.output_dir / "final_AMP_DE_results.tsv", final_rows)

    summary = []
    for contrast in contrasts:
        for family in ("ALL", "Defensin", "Snakin_GASA", "nsLTP"):
            selected = [
                row for row in final_rows
                if row["contrast"] == contrast and (family == "ALL" or row["family"] == family)
            ]
            quantified = [row for row in selected if row["quantification_status"] == "QUANTIFIED"]
            significant = [row for row in quantified if de_call(str(row["DE_call_FDR_lt_0.05_abs_log2FC_gt_1"]))]
            directions = Counter(
                "up" if float(row["log2FoldChange"]) > 1 else "down"
                for row in significant
            )
            summary.append(
                {
                    "contrast": contrast,
                    "family": family,
                    "primary_catalogue_members": len(selected),
                    "quantified_members": len(quantified),
                    "retained_for_DE_members": sum(row["prefilter_status"] == "RETAINED_TOTAL_COUNT_GE_10" for row in selected),
                    "filtered_all_zero_members": sum(row["prefilter_status"] == "FILTERED_ALL_ZERO_COUNTS" for row in selected),
                    "filtered_low_total_count_members": sum(row["prefilter_status"] == "FILTERED_TOTAL_COUNT_LT_10" for row in selected),
                    "members_with_adjusted_P": sum(row["DE_test_status"] == "TESTED_WITH_ADJUSTED_P" for row in selected),
                    "DE_members": len(significant),
                    "upregulated": directions["up"],
                    "downregulated": directions["down"],
                    "threshold": "Benjamini-Hochberg FDR < 0.05 and |unshrunken log2FC| > 1",
                }
            )
    write_tsv(args.output_dir / "final_AMP_DE_summary.tsv", summary)

    candidates = [
        {
            **row,
            "candidate_evidence_label": "salt-responsive transcriptome-supported candidate",
            "functional_validation_status": "NOT_EXPERIMENTALLY_VALIDATED",
        }
        for row in final_rows
        if row["contrast"] == "S400_vs_CK"
        and row["quantification_status"] == "QUANTIFIED"
        and de_call(str(row["DE_call_FDR_lt_0.05_abs_log2FC_gt_1"]))
    ]
    candidates.sort(key=lambda row: (float(row["padj_BH"]), -abs(float(row["log2FoldChange"]))))
    if candidates:
        write_tsv(args.output_dir / "S400_candidate_priority.tsv", candidates)
    else:
        (args.output_dir / "S400_candidate_priority.tsv").write_text(
            "status\treason\nNO_CANDIDATES\tNo primary-catalogue member met the predeclared S400-versus-CK threshold\n",
            encoding="utf-8",
        )

    matrix_subset(args.tpm, catalogue, args.output_dir / "final_AMP_TPM.tsv", audit, False)
    matrix_subset(args.normalized_counts, catalogue, args.output_dir / "final_AMP_normalized_counts.tsv", audit, True)
    matrix_subset(args.vst, catalogue, args.output_dir / "final_AMP_VST.tsv", audit, True)
    (args.output_dir / "FINAL_EXPRESSION_CATALOGUE_COMPLETE.PASS").write_text(
        f"PASS\tprimary_catalogue={len(catalogue)}\tS400_DE={len(candidates)}\n",
        encoding="utf-8",
    )
    print(f"PASS primary catalogue={len(catalogue)} S400 DE={len(candidates)}")


if __name__ == "__main__":
    main()
