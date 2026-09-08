#!/usr/bin/env python3
"""Prepare and summarize a genome-level Tamarix defensin homology rescue."""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path


ORIGINAL = Path("/path/to/rsoamp")
REVISION = ORIGINAL / "revision_R1_20260902"
OUT = REVISION / "05_comparative_genomics"
PRED = OUT / "predictions"
GFF = OUT / "inputs/Tamarix_austromongolica/tau.longest.gff3"


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
    return {key: "".join(value).replace("*", "") for key, value in records.items()}


def prepare() -> None:
    with (ORIGINAL / "tables/T1_AMP_members.csv").open(encoding="utf-8-sig") as handle:
        members = [
            row for row in csv.DictReader(handle)
            if row["family"] == "Defensin" and row["source"] == "annotated"
        ]
    sequences = read_fasta(ORIGINAL / "intermediate/members_final.faa")
    with (PRED / "rso_annotated_defensins.faa").open("w", encoding="utf-8") as handle:
        for member in members:
            member_id = member["member_id"]
            handle.write(f">RSO_{member_id}\n{sequences[member_id]}\n")
    print(f"Prepared {len(members)} annotated R. soongarica defensin queries")


def parse_attributes(text: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in text.split(";") if "=" in part)


def load_genes() -> dict[str, list[tuple[int, int, str]]]:
    genes = defaultdict(list)
    with GFF.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9 and fields[2] == "gene":
                attrs = parse_attributes(fields[8])
                genes[fields[0]].append((int(fields[3]), int(fields[4]), attrs.get("ID", "")))
    return genes


def parse() -> None:
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "qlen",
    ]
    hsps = []
    with (PRED / "tamarix_defensin_tblastn.tsv").open(encoding="utf-8") as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if not fields:
                continue
            row = dict(zip(columns, fields))
            if int(row["length"]) >= 20:
                hsps.append(row)

    grouped = defaultdict(list)
    for row in hsps:
        strand = "+" if int(row["sstart"]) <= int(row["send"]) else "-"
        start, end = sorted((int(row["sstart"]), int(row["send"])))
        grouped[(row["sseqid"], strand)].append((start, end, row))

    clusters = []
    for (chrom, strand), items in grouped.items():
        items.sort(key=lambda item: (item[0], item[1], item[2]["qseqid"]))
        current = []
        current_end = -1
        for start, end, row in items:
            if current and start - current_end > 10_000:
                clusters.append((chrom, strand, current))
                current = []
                current_end = -1
            current.append((start, end, row))
            current_end = max(current_end, end) if current_end >= 0 else end
        if current:
            clusters.append((chrom, strand, current))

    genes = load_genes()
    output = []
    for index, (chrom, strand, items) in enumerate(sorted(clusters), start=1):
        start = min(item[0] for item in items)
        end = max(item[1] for item in items)
        overlaps = [gene_id for gene_start, gene_end, gene_id in genes.get(chrom, []) if gene_start <= end and gene_end >= start]
        queries = sorted({item[2]["qseqid"] for item in items})
        output.append(
            {
                "cluster_id": f"TAU_DEF_TBLASTN_{index:03d}",
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "span_bp": end - start + 1,
                "n_hsps": len(items),
                "n_queries": len(queries),
                "query_ids": ";".join(queries),
                "best_evalue": min(float(item[2]["evalue"]) for item in items),
                "best_bitscore": max(float(item[2]["bitscore"]) for item in items),
                "max_query_coverage": max(
                    (int(item[2]["qend"]) - int(item[2]["qstart"]) + 1) / int(item[2]["qlen"])
                    for item in items
                ),
                "overlapping_annotated_genes": ";".join(sorted(set(overlaps))),
                "annotation_status": "ANNOTATED_GENE_OVERLAP" if overlaps else "INTERGENIC_HOMOLOGY_CANDIDATE",
                "interpretation": "Homology locus only; not counted as a gene without an intact model and transcript evidence.",
            }
        )

    target = OUT / "tamarix_defensin_genome_rescue.tsv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(
        f"Parsed {len(hsps)} HSPs into {len(output)} loci; "
        f"{sum(row['annotation_status'] == 'INTERGENIC_HOMOLOGY_CANDIDATE' for row in output)} are intergenic"
    )

    translated = PRED / "tamarix_defensin_miniprot_translated.txt"
    accepted = []
    if translated.exists():
        pending = None
        for line in translated.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            if line.startswith("##STA\t") and pending is not None:
                sequence = line.split("\t", 1)[1].rstrip("*")
                query_coverage = (pending[3] - pending[2]) / pending[1]
                tags = set(pending[12:])
                if query_coverage >= 0.9 and "fs:i:0" in tags and "st:i:0" in tags and sequence.startswith("M"):
                    accepted.append(
                        {
                            "query_id": pending[0],
                            "query_coverage": query_coverage,
                            "chrom": pending[5],
                            "start": int(pending[7]) + 1,
                            "end": int(pending[8]),
                            "strand": pending[4],
                            "alignment_score": next(tag.split(":")[-1] for tag in pending[12:] if tag.startswith("AS:i:")),
                            "sequence": sequence,
                        }
                    )
                pending = None
            elif not line.startswith("#"):
                fields = line.split("\t")
                pending = [fields[0], int(fields[1]), int(fields[2]), int(fields[3]), *fields[4:]]

    unique = []
    seen = set()
    for row in accepted:
        if row["sequence"] in seen:
            continue
        seen.add(row["sequence"])
        unique.append(row)
    for index, row in enumerate(unique, start=1):
        row["candidate_id"] = f"TAU_DEF_MP{index:03d}"

    fasta = PRED / "tamarix_Defensin_high_confidence.faa"
    with fasta.open("w", encoding="utf-8") as handle:
        for row in unique:
            handle.write(f">TAU|{row['candidate_id']}|Defensin\n{row['sequence']}\n")
    manifest = OUT / "tamarix_defensin_miniprot_models.tsv"
    if unique:
        columns = [key for key in unique[0] if key != "sequence"] + ["sequence"]
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(unique)
    print(f"Retained {len(unique)} unique complete miniprot defensin homology models for phylogenetic context")


def main() -> int:
    mode = sys.argv[1]
    if mode == "prepare":
        prepare()
    elif mode == "parse":
        parse()
    else:
        raise ValueError(mode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
