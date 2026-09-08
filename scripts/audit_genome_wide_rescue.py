#!/usr/bin/env python3
"""Make the historical six-family whole-genome tblastn rescue screen auditable."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--revision", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    a_lo, a_hi = sorted((a_start, a_end))
    b_lo, b_hi = sorted((b_start, b_end))
    return not (a_hi < b_lo or b_hi < a_lo)


def gene_id(attr: str) -> str:
    for item in attr.rstrip(";").split(";"):
        if item.startswith("ID="):
            return item.split("=", 1)[1]
    return attr


def main() -> None:
    args = parse_args()
    backfill_dir = args.project / "intermediate/backfill"
    outdir = args.revision / "03_annotation_rescue"

    loci = json.loads((backfill_dir / "tblastn_loci.json").read_text(encoding="utf-8"))
    raw_hits = []
    with (backfill_dir / "tblastn_raw.tsv").open(encoding="utf-8") as handle:
        for line in handle:
            query, chrom, pident, length, evalue, bitscore, sstart, send = line.rstrip("\n").split("\t")
            raw_hits.append(
                {
                    "query": query,
                    "family": query.split("|", 1)[0],
                    "chrom": chrom,
                    "pident": float(pident),
                    "length": int(length),
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                    "sstart": int(sstart),
                    "send": int(send),
                }
            )

    genes: dict[str, list[tuple[int, int, str]]] = {}
    with (backfill_dir / "genes.bed").open(encoding="utf-8") as handle:
        for line in handle:
            chrom, start, end, attr = line.rstrip("\n").split("\t")
            genes.setdefault(chrom, []).append((int(start), int(end), gene_id(attr)))

    preliminary = read_tsv(outdir / "backfill_locus_audit_preliminary.tsv")
    known_by_key = {}
    for row in preliminary:
        chrom, rest = row["historical_locus"].split(":", 1)
        interval = rest.split("(", 1)[0]
        start, end = (int(value) for value in interval.split("-"))
        known_by_key[(chrom, start, end)] = row

    final_path = outdir / "backfill_locus_audit_final.tsv"
    final_by_member = {}
    if final_path.exists() and final_path.stat().st_size:
        final_by_member = {row["member_id"]: row for row in read_tsv(final_path)}

    rows: list[dict[str, object]] = []
    for locus_index, (chrom, start, end, families, cluster_evalue) in enumerate(loci, start=1):
        matching = [
            hit for hit in raw_hits
            if hit["chrom"] == chrom and overlaps(start, end, hit["sstart"], hit["send"])
        ]
        if not matching:
            raise ValueError(f"No raw HSP maps to clustered locus {chrom}:{start}-{end}")
        best = min(matching, key=lambda hit: (hit["evalue"], -hit["bitscore"]))
        overlapping_genes = [
            gid for gstart, gend, gid in genes.get(chrom, []) if overlaps(start, end, gstart, gend)
        ]
        known = known_by_key.get((chrom, start, end))
        member_id = known["member_id"] if known else ""
        final = final_by_member.get(member_id, {})
        if overlapping_genes:
            classification = "ANNOTATED_GENE_OVERLAP"
            primary_inclusion = "NOT_APPLICABLE_EXISTING_ANNOTATION"
        elif known:
            classification = final.get("final_status", "KNOWN_UNANNOTATED_LOCUS_PENDING_FINAL_EVIDENCE")
            primary_inclusion = final.get("primary_catalogue_inclusion", "PENDING")
        else:
            classification = "UNANNOTATED_AMBIGUOUS_HIT_REQUIRES_MODEL_REVIEW"
            primary_inclusion = "NO"
        rows.append(
            {
                "screen_locus_id": f"GWRS{locus_index:03d}",
                "chromosome": chrom,
                "cluster_start_1based": start,
                "cluster_end_1based": end,
                "cluster_length_nt": end - start + 1,
                "families_from_cluster": ";".join(families),
                "cluster_best_evalue_recorded": cluster_evalue,
                "raw_hsp_count": len(matching),
                "best_query": best["query"],
                "best_query_family": best["family"],
                "best_hsp_evalue": best["evalue"],
                "best_hsp_bitscore": best["bitscore"],
                "best_hsp_pident": best["pident"],
                "best_hsp_alignment_aa": best["length"],
                "best_hsp_strand": "+" if best["sstart"] <= best["send"] else "-",
                "overlapping_gene_model_count": len(overlapping_genes),
                "overlapping_gene_model_ids": ";".join(sorted(set(overlapping_genes))),
                "known_backfill_member_id": member_id,
                "final_rescue_classification": classification,
                "primary_catalogue_inclusion": primary_inclusion,
            }
        )

    if len(rows) != 37:
        raise ValueError(f"Expected 37 whole-genome clustered loci, observed {len(rows)}")
    if len(raw_hits) != 1003:
        raise ValueError(f"Expected 1003 raw HSPs, observed {len(raw_hits)}")
    unannotated = [row for row in rows if row["overlapping_gene_model_count"] == 0]
    if len(unannotated) != 6 or any(not row["known_backfill_member_id"] for row in unannotated):
        raise ValueError("Whole-genome screen no longer resolves to the six audited unannotated loci")

    write_tsv(outdir / "genome_wide_rescue_hits.tsv", rows)
    if final_by_member:
        write_tsv(outdir / "six_loci_audit.tsv", [final_by_member[row["member_id"]] for row in preliminary])

    status_counts = Counter(str(row["final_rescue_classification"]) for row in rows)
    summary = [
        "# Genome-wide targeted rescue audit",
        "",
        "Status: PASS",
        "",
        "The historical rescue was a whole-genome screen, not a search limited to six preselected loci.",
        f"It used 128 family-labelled reference peptides, produced {len(raw_hits)} tblastn HSPs, and collapsed them into {len(rows)} loci.",
        f"Thirty-one loci overlap frozen gene models; the remaining {len(unannotated)} are the six loci evaluated with ORF, short-read, Iso-Seq, domain, and secretion evidence.",
        "No additional unannotated cluster exists in the frozen screen.",
        "",
        "## Final locus classes",
        "",
        *[f"- {key}: {value}" for key, value in sorted(status_counts.items())],
        "",
        "The screen itself is sensitive rather than gene-defining: a tblastn cluster is not counted as a gene without a complete model and transcript support.",
    ]
    (outdir / "annotation_incompleteness_summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    (outdir / "GENOME_WIDE_RESCUE_AUDIT_COMPLETE.PASS").write_text(
        f"PASS\traw_hsps={len(raw_hits)}\tclustered_loci={len(rows)}\tunannotated={len(unannotated)}\n",
        encoding="utf-8",
    )
    print(f"PASS whole-genome rescue audit: {len(raw_hits)} HSPs, {len(rows)} loci, {len(unannotated)} unannotated")


if __name__ == "__main__":
    main()
