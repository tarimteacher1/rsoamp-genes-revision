#!/usr/bin/env python3
"""Build the primary AMP gene catalogue after nsLTP and backfill audits."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--revision", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_attrs(text: str) -> dict[str, str]:
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def transcript_to_gene(gff: Path) -> dict[str, str]:
    mapping = {}
    with gff.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("ID") and attrs.get("Parent"):
                mapping[attrs["ID"]] = attrs["Parent"].split(",", 1)[0]
    return mapping


def read_fasta(path: Path) -> dict[str, str]:
    sequences = {}
    current = None
    chunks = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        sequences[current] = "".join(chunks)
    return sequences


def main() -> None:
    args = parse_args()
    outdir = args.revision / "03_annotation_rescue"
    nsdir = args.revision / "01_nslTP_curation"
    unmapped = args.revision / "04_unmapped_reads"

    t1 = read_table(args.project / "tables/T1_AMP_members.csv", ",")
    ns_rows = read_table(nsdir / "nsLTP_evidence_matrix.tsv")
    ns_by_member = {row["member_id"]: row for row in ns_rows if row["member_id"]}
    backfill = read_table(outdir / "backfill_locus_audit_preliminary.tsv")
    short_support = {row["member_id"]: row for row in read_table(unmapped / "recovered_loci_read_support.tsv")}
    long_support = {row["member_id"]: row for row in read_table(outdir / "isoseq_backfill_support.tsv")}
    isoseq_models = {row["member_id"]: row for row in read_table(outdir / "isoseq_backfill_gene_models.tsv")}
    tx2gene = transcript_to_gene(args.project / "data/genome.gff")
    proteins = read_fasta(args.project / "data/protein.fa")

    if len(t1) != 73 or len(backfill) != 6:
        raise ValueError(f"Unexpected historical input sizes: T1={len(t1)}, backfill={len(backfill)}")
    if (
        set(row["member_id"] for row in backfill) != set(short_support)
        or set(short_support) != set(long_support)
        or set(long_support) != set(isoseq_models)
    ):
        raise ValueError("Backfill support tables do not contain the same six member IDs")

    backfill_final = []
    for row in backfill:
        member = row["member_id"]
        short = short_support[member]
        long = long_support[member]
        isoseq_model = isoseq_models[member]
        local_complete_orf = row["upstream_inframe_methionine"] == "YES" and row["downstream_stop_before_next_methionine"] == "YES"
        isoseq_complete_model = isoseq_model["complete_isoseq_supported_gene_model_gate"] == "PASS"
        complete_orf_gate = local_complete_orf or isoseq_complete_model
        annotation_gate = bool(row["overlapping_gff_ids"]) or isoseq_complete_model
        transcript_gate = (
            short["short_read_evidence_class"] == "CONSISTENT_SHORT_READ_LOCUS_SUPPORT"
            and long["transcript_evidence_class"] == "STRONG_CONTIGUOUS_TRANSCRIPT_SUPPORT"
            and isoseq_complete_model
        )
        include = complete_orf_gate and annotation_gate and transcript_gate
        final_status = "CONFIRMED_COMPLETE_GENE_MODEL" if include else "GENOMIC_AMP_LIKE_ORF_CANDIDATE_NOT_CONFIRMED_AS_GENE"
        reason_parts = []
        if not complete_orf_gate:
            reason_parts.append("neither the local genomic interval nor an Iso-Seq transcript model supplied a complete start-to-stop ORF")
        if not annotation_gate:
            reason_parts.append("no overlapping frozen GFF model or complete Iso-Seq rescue model")
        if not transcript_gate:
            reason_parts.append("combined short- and long-read gates for a complete transcript were not met")
        backfill_final.append(
            {
                **row,
                "short_read_evidence_class": short["short_read_evidence_class"],
                "short_read_pooled_breadth_depth_ge1": short["pooled_breadth_depth_ge1"],
                "short_read_supporting_samples": short["samples_with_fragment_and_ge50pct_breadth_support"],
                "isoseq_transcript_evidence_class": long["transcript_evidence_class"],
                "isoseq_maximum_aligned_coding_fraction": long["maximum_aligned_coding_fraction"],
                "isoseq_complete_model_transcript_id": isoseq_model["best_transcript_id"],
                "isoseq_complete_model_coordinates": isoseq_model["best_transcript_coordinates"],
                "isoseq_complete_orf_length_aa": isoseq_model["best_complete_orf_length_aa"],
                "isoseq_complete_orf_protein_sequence": isoseq_model["best_complete_orf_protein_sequence"],
                "complete_orf_gate": "PASS" if complete_orf_gate else "FAIL",
                "annotation_model_gate": "PASS" if annotation_gate else "FAIL",
                "combined_transcript_support_gate": "PASS" if transcript_gate else "FAIL",
                "final_gene_model_status": final_status,
                "primary_catalogue_inclusion": "YES" if include else "NO",
                "exclusion_reason": "; ".join(reason_parts),
            }
        )
    write_tsv(outdir / "backfill_locus_audit_final.tsv", backfill_final)

    backfill_by_member = {row["member_id"]: row for row in backfill_final}
    catalogue = []
    for row in t1:
        member = row["member_id"]
        if row["family"] == "nsLTP":
            evidence = ns_by_member[member]
            included = evidence["final_classification"] == "A_CANONICAL_NSLTP"
            status = "INCLUDED_ANNOTATED_GENE" if included else "EXCLUDED_AMBIGUOUS_NSLTP_PROLAMIN_BOUNDARY"
            final_subclass = "canonical_nsLTP" if included else "ambiguous_nsLTP_prolamin_boundary"
            final_subtype = evidence["final_subtype"]
            confidence = evidence["classification_confidence"]
            reason = evidence["classification_rule_trace"]
        elif row["source"] == "tblastn_backfill":
            audit = backfill_by_member[member]
            included = audit["primary_catalogue_inclusion"] == "YES"
            status = "INCLUDED_RESCUED_COMPLETE_GENE" if included else "EXCLUDED_GENOMIC_ORF_CANDIDATE"
            final_subclass = "confirmed_rescued_gene" if included else "genomic_AMP_like_ORF_candidate"
            final_subtype = row["subtype"]
            confidence = "HIGH" if included else "HIGH_FOR_EXCLUSION_FROM_GENE_COUNT"
            reason = audit["exclusion_reason"] if not included else "All complete-ORF, annotation, and transcript-support gates passed"
        else:
            included = True
            status = "INCLUDED_ANNOTATED_GENE"
            final_subclass = "canonical_family_member"
            final_subtype = row["subtype"]
            confidence = "HIGH"
            reason = "Annotated gene model with family-specific domain evidence in the audited original workflow"

        protein_id = row["protein_id"]
        if row["source"] == "annotated":
            gene_id = tx2gene.get(protein_id, "")
        elif included:
            gene_id = f"ISOSEQ_RESCUE_{member}"
            protein_id = f"ISOSEQ_RESCUE_{member}.p1"
            proteins[protein_id] = backfill_by_member[member]["isoseq_complete_orf_protein_sequence"]
        else:
            gene_id = ""
        if included and row["source"] == "annotated" and (not gene_id or protein_id not in proteins):
            raise ValueError(f"Included annotated member lacks GFF/protein mapping: {member} {protein_id}")
        catalogue.append(
            {
                "member_id": member,
                "family": row["family"],
                "protein_id": protein_id,
                "gene_id": gene_id,
                "chrom": row["chrom"],
                "start": row["start"],
                "end": row["end"],
                "strand": row["strand"],
                "length_aa": row["length_aa"],
                "historical_source": row["source"],
                "historical_subclass": row["subclass"],
                "final_subclass": final_subclass,
                "final_subtype": final_subtype,
                "catalogue_status": status,
                "classification_confidence": confidence,
                "classification_reason": reason,
            }
        )

    write_tsv(outdir / "final_amp_catalogue_all_audited_members.tsv", catalogue)
    included_catalogue = [row for row in catalogue if row["catalogue_status"].startswith("INCLUDED")]
    write_tsv(outdir / "final_amp_primary_catalogue.tsv", included_catalogue)

    counts = Counter(row["family"] for row in included_catalogue)
    count_rows = [
        {"family": family, "primary_gene_count": counts.get(family, 0), "count_definition": "Included complete annotated or fully rescued gene models"}
        for family in ("Defensin", "Snakin_GASA", "nsLTP")
    ]
    count_rows.append(
        {"family": "TOTAL", "primary_gene_count": len(included_catalogue), "count_definition": "Sum of the three primary AMP gene families"}
    )
    write_tsv(outdir / "final_amp_family_counts.tsv", count_rows)

    with (outdir / "final_amp_primary_catalogue.faa").open("w", encoding="ascii") as handle:
        for row in included_catalogue:
            sequence = proteins[row["protein_id"]]
            handle.write(f">{row['member_id']} protein_id={row['protein_id']} gene_id={row['gene_id']} family={row['family']}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    report = [
        "# Final AMP catalogue completeness report",
        "",
        "Status: PASS",
        "",
        f"The primary gene catalogue contains {len(included_catalogue)} complete annotated or fully rescued gene models: "
        + ", ".join(f"{family}={counts.get(family, 0)}" for family in ("Defensin", "Snakin_GASA", "nsLTP"))
        + ".",
        f"All {sum(row['historical_source'] == 'tblastn_backfill' for row in catalogue)} historical tblastn backfills were assessed against exact genomic encoding, local start/stop context, frozen GFF overlap, nine-sample short-read support, and public Iso-Seq support.",
        f"Backfill candidates retained in the primary gene count: {sum(row['primary_catalogue_inclusion'] == 'YES' for row in backfill_final)}.",
        "Excluded backfill sequences remain reported as genomic AMP-like ORF candidates and are not silently discarded, but they are not used for gene-family counts, expression, synteny, promoter, or duplication analyses.",
        "The catalogue is evidence-defined; historical totals were not used as retention targets.",
    ]
    (outdir / "catalogue_completeness_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (outdir / "FINAL_AMP_CATALOGUE_COMPLETE.PASS").write_text(
        f"PASS\tprimary_members={len(included_catalogue)}\tbackfill_included={sum(row['primary_catalogue_inclusion'] == 'YES' for row in backfill_final)}\n",
        encoding="utf-8",
    )
    print("\n".join(report))


if __name__ == "__main__":
    main()
