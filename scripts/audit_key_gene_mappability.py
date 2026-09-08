#!/usr/bin/env python3
"""Audit AMP transcript similarity and Salmon assignment ambiguity."""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--quant-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--makeblastdb", default="/usr/bin/makeblastdb")
    parser.add_argument("--blastn", default="/usr/bin/blastn")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
    return {key: "".join(value) for key, value in records.items()}


def salmon_ambiguity(quant_root: Path) -> dict[str, list[float]]:
    values: dict[str, list[float]] = {}
    for run_dir in sorted(path for path in quant_root.iterdir() if path.is_dir()):
        quant_path = run_dir / "quant.sf"
        ambiguity_path = run_dir / "aux_info/ambig_info.tsv"
        if not quant_path.exists() or not ambiguity_path.exists():
            continue
        quant = read_tsv(quant_path)
        ambiguity = read_tsv(ambiguity_path)
        if len(quant) != len(ambiguity):
            raise ValueError(f"Salmon quant/ambiguity row mismatch for {run_dir.name}")
        for quant_row, ambiguity_row in zip(quant, ambiguity):
            transcript = ambiguity_row.get("Name") or ambiguity_row.get("Transcript") or quant_row["Name"]
            unique = float(ambiguity_row["UniqueCount"])
            ambiguous = float(ambiguity_row["AmbigCount"])
            denominator = unique + ambiguous
            values.setdefault(transcript, []).append(ambiguous / denominator if denominator else 0.0)
    return values


def main() -> None:
    args = parse_args()
    catalogue = read_tsv(args.catalogue)
    transcripts = read_fasta(args.transcripts)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    query_path = args.output.parent / "primary_AMP_transcripts.fasta"
    missing = []
    with query_path.open("w", encoding="ascii") as handle:
        for row in catalogue:
            transcript = row["protein_id"]
            if transcript not in transcripts:
                missing.append(transcript)
                continue
            handle.write(f">{transcript}\n{transcripts[transcript]}\n")
    if not query_path.stat().st_size:
        raise ValueError("No primary AMP transcript was found in the reconstructed transcriptome")

    database = args.output.parent / "full_transcriptome_blastdb"
    subprocess.run(
        [args.makeblastdb, "-in", str(args.transcripts), "-dbtype", "nucl", "-parse_seqids", "-out", str(database)],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    blast_path = args.output.parent / "primary_AMP_vs_full_transcriptome.blastn.tsv"
    with blast_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                args.blastn,
                "-task", "blastn",
                "-query", str(query_path),
                "-db", str(database),
                "-evalue", "1e-10",
                "-dust", "yes",
                "-max_target_seqs", "500",
                "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
            ],
            check=True,
            stdout=handle,
        )

    hits: dict[str, list[dict[str, object]]] = {}
    with blast_path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore in reader:
            if qseqid == sseqid:
                continue
            alignment = int(length)
            query_length = int(qlen)
            hits.setdefault(qseqid, []).append(
                {
                    "subject": sseqid,
                    "identity": float(pident),
                    "aligned_length": alignment,
                    "query_coverage": alignment / query_length,
                    "bitscore": float(bitscore),
                }
            )

    ambiguity = salmon_ambiguity(args.quant_root)
    rows = []
    for member in catalogue:
        transcript = member["protein_id"]
        sequence = transcripts.get(transcript)
        transcript_hits = sorted(hits.get(transcript, []), key=lambda row: row["bitscore"], reverse=True)
        best = transcript_hits[0] if transcript_hits else None
        close_hits = [row for row in transcript_hits if row["identity"] >= 95 and row["query_coverage"] >= 0.8]
        fractions = ambiguity.get(transcript, [])
        median_ambiguity = statistics.median(fractions) if fractions else None
        sequence_risk = bool(close_hits)
        assignment_risk = median_ambiguity is not None and median_ambiguity >= 0.5
        if sequence is None:
            risk = "NOT_ASSESSED_CUSTOM_RESCUE_TRANSCRIPT"
        elif sequence_risk and assignment_risk:
            risk = "HIGH"
        elif sequence_risk or assignment_risk or (median_ambiguity is not None and median_ambiguity >= 0.2):
            risk = "MODERATE"
        else:
            risk = "LOW"
        rows.append(
            {
                "member_id": member["member_id"],
                "family": member["family"],
                "gene_id": member["gene_id"],
                "transcript_id": transcript,
                "transcript_length_nt": len(sequence) if sequence else "",
                "best_nonself_transcript": best["subject"] if best else "",
                "best_nonself_identity_pct": f"{best['identity']:.3f}" if best else "",
                "best_nonself_query_coverage": f"{best['query_coverage']:.4f}" if best else "",
                "near_identical_nonself_hits_95pct_80pct_coverage": len(close_hits),
                "samples_with_salmon_ambiguity_metrics": len(fractions),
                "median_salmon_ambiguous_fraction": f"{median_ambiguity:.4f}" if median_ambiguity is not None else "",
                "mappability_risk": risk,
                "interpretation": "Sequence-similarity and read-assignment audit; not a qPCR primer-specificity test",
            }
        )

    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    missing_path = args.output.parent / "primary_AMP_transcripts_missing_from_frozen_annotation.tsv"
    missing_path.write_text(
        "transcript_id\treason\n"
        + "".join(f"{transcript}\tcustom rescue transcript absent from frozen GFF-derived Salmon index\n" for transcript in missing),
        encoding="utf-8",
    )
    print(f"PASS mappability audit members={len(rows)} missing_transcripts={len(missing)}")


if __name__ == "__main__":
    main()
