#!/usr/bin/env python3
"""Test whether target Iso-Seq transcript models contain complete historical AMP ORFs."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
BACKFILL_IDS = {"RsDEF2", "RsDEF3", "RsDEF4", "RsDEF5", "RsDEF7", "RsSNA13"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-peptides", required=True, type=Path)
    parser.add_argument("--backfill-audit", required=True, type=Path)
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    records = {}
    current = None
    chunks = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks).upper().replace("*", "")
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks).upper().replace("*", "")
    return records


def translate(sequence: str, frame: int) -> str:
    return "".join(CODON_TABLE.get(sequence[index : index + 3], "X") for index in range(frame, len(sequence) - 2, 3))


def parse_gtf(path: Path) -> dict[str, dict[str, object]]:
    models: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not match:
                continue
            transcript_id = match.group(1)
            models[transcript_id] = {
                "chrom": fields[0], "start": int(fields[3]), "end": int(fields[4]), "strand": fields[6]
            }
    return models


def complete_orfs_containing(sequence: str, peptide: str) -> list[dict[str, object]]:
    results = []
    for frame in range(3):
        translation = translate(sequence, frame)
        peptide_start = translation.find(peptide)
        while peptide_start >= 0:
            prior_stop = translation.rfind("*", 0, peptide_start)
            start_m = translation.rfind("M", prior_stop + 1, peptide_start + 1)
            next_stop = translation.find("*", peptide_start + len(peptide))
            complete = start_m >= 0 and next_stop >= 0
            complete_protein = translation[start_m:next_stop] if complete else ""
            results.append(
                {
                    "frame": frame,
                    "peptide_aa_start": peptide_start + 1,
                    "orf_start_aa": start_m + 1 if start_m >= 0 else "",
                    "orf_stop_aa": next_stop + 1 if next_stop >= 0 else "",
                    "complete_start_to_stop_orf": "YES" if complete else "NO",
                    "complete_orf_length_aa": next_stop - start_m if complete else "",
                    "complete_orf_protein_sequence": complete_protein,
                }
            )
            peptide_start = translation.find(peptide, peptide_start + 1)
    return results


def main() -> None:
    args = parse_args()
    peptides = {key: value for key, value in read_fasta(args.historical_peptides).items() if key in BACKFILL_IDS}
    transcripts = read_fasta(args.transcripts)
    models = parse_gtf(args.gtf)
    with args.backfill_audit.open(encoding="utf-8", newline="") as handle:
        loci = {row["member_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(peptides) != BACKFILL_IDS or set(loci) != BACKFILL_IDS:
        raise ValueError("Historical peptide or locus set does not contain the six expected members")

    output_rows = []
    for member_id in sorted(BACKFILL_IDS):
        peptide = peptides[member_id]
        matches = []
        for transcript_id, sequence in transcripts.items():
            for orf in complete_orfs_containing(sequence, peptide):
                matches.append((transcript_id, sequence, orf))
        complete_matches = [item for item in matches if item[2]["complete_start_to_stop_orf"] == "YES"]
        best_pool = complete_matches or matches
        best = max(best_pool, key=lambda item: int(item[2]["complete_orf_length_aa"] or 0), default=None)
        if best:
            transcript_id, sequence, orf = best
            model = models.get(transcript_id, {})
            transcript_coordinates = (
                f"{model.get('chrom')}:{model.get('start')}-{model.get('end')}({model.get('strand')})" if model else ""
            )
        else:
            transcript_id, sequence, orf, transcript_coordinates = "", "", {}, ""
        output_rows.append(
            {
                "member_id": member_id,
                "family": loci[member_id]["family"],
                "historical_locus": loci[member_id]["historical_locus"],
                "historical_peptide_length_aa": len(peptide),
                "stringtie_transcripts_containing_exact_peptide": len({item[0] for item in matches}),
                "stringtie_transcripts_with_complete_start_stop_orf_containing_peptide": len({item[0] for item in complete_matches}),
                "best_transcript_id": transcript_id,
                "best_transcript_coordinates": transcript_coordinates,
                "best_transcript_length_nt": len(sequence) if sequence else "",
                "best_orf_frame": orf.get("frame", ""),
                "best_orf_start_aa": orf.get("orf_start_aa", ""),
                "best_orf_stop_aa": orf.get("orf_stop_aa", ""),
                "best_complete_orf_length_aa": orf.get("complete_orf_length_aa", ""),
                "best_complete_orf_protein_sequence": orf.get("complete_orf_protein_sequence", ""),
                "complete_isoseq_supported_gene_model_gate": "PASS" if complete_matches else "FAIL",
                "interpretation_limit": "StringTie reconstruction from target-filtered public Iso-Seq alignments requires domain and locus-consistency review before gene acceptance",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {len(output_rows)} Iso-Seq gene-model audit rows")


if __name__ == "__main__":
    main()
