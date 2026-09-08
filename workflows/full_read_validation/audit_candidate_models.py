#!/usr/bin/env python3
"""Verify candidate ORFs, splice motifs, and genomic CDS coordinates."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')
STOPS = {"TAA", "TAG", "TGA"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                current = identifier if identifier in wanted else None
                if current is not None:
                    records[current] = []
            elif current is not None:
                records[current].append(line)
    return {key: "".join(value).upper() for key, value in records.items()}


def parse_gtf(path: Path) -> dict[str, dict[str, object]]:
    models: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            attrs = dict(ATTRIBUTE_RE.findall(fields[8]))
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            model = models.setdefault(tid, {"chrom": fields[0], "strand": fields[6], "exons": []})
            if fields[2] == "exon":
                model["exons"].append((int(fields[3]) - 1, int(fields[4])))
    for model in models.values():
        model["exons"] = sorted(set(model["exons"]))
    return models


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(sequence: str) -> str:
    return "".join(CODONS.get(sequence[index : index + 3], "X") for index in range(0, len(sequence), 3))


def project_interval(exons: list[tuple[int, int]], strand: str, start0: int, end0: int) -> list[tuple[int, int]]:
    ordered = exons if strand == "+" else list(reversed(exons))
    projected: list[tuple[int, int]] = []
    offset = 0
    for exon_start, exon_end in ordered:
        length = exon_end - exon_start
        local_start = max(0, start0 - offset)
        local_end = min(length, end0 - offset)
        if local_start < local_end:
            if strand == "+":
                projected.append((exon_start + local_start, exon_start + local_end))
            else:
                projected.append((exon_end - local_end, exon_end - local_start))
        offset += length
    if sum(end - start for start, end in projected) != end0 - start0:
        raise ValueError("Spliced interval projection did not preserve length")
    return projected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-evidence", required=True, type=Path)
    parser.add_argument("--candidate-gtf", required=True, type=Path)
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--priority-peptides", required=True, type=Path)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    candidates = [
        row for row in read_tsv(args.candidate_evidence)
        if row["model_origin"] == "INTERGENIC_CANDIDATE" and row["review_priority"] == "FULL_READ_MODEL_VALIDATION_PRIORITY"
    ]
    tids = {row["transcript_id"] for row in candidates}
    transcripts = read_fasta(args.transcripts, tids)
    peptides = {
        key.removesuffix("|start_stop"): value
        for key, value in read_fasta(args.priority_peptides, {row["orf_id"] + "|start_stop" for row in candidates}).items()
    }
    models = parse_gtf(args.candidate_gtf)
    if set(transcripts) != tids or set(peptides) != {row["orf_id"] for row in candidates}:
        raise ValueError("Candidate sequence identities do not join")

    import pysam

    genome = pysam.FastaFile(str(args.genome))
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    audit: list[dict[str, object]] = []
    gff_lines = ["##gff-version 3\n"]
    for index, row in enumerate(sorted(candidates, key=lambda value: (value["chrom"], int(value["start"]))), start=1):
        tid, oid = row["transcript_id"], row["orf_id"]
        model, transcript, peptide = models[tid], transcripts[tid], peptides[oid]
        first_m = int(row["first_methionine_position_aa"])
        cds_start0 = int(row["orf_start_in_transcript"]) - 1 + (first_m - 1) * 3
        cds_end0 = cds_start0 + len(peptide) * 3
        coding_nt = transcript[cds_start0:cds_end0]
        stop = transcript[cds_end0 : cds_end0 + 3]
        translated = translate(coding_nt)
        if translated != peptide or stop not in STOPS:
            raise ValueError(f"ORF projection/translation mismatch: {oid}")
        exons = model["exons"]
        projected = project_interval(exons, str(model["strand"]), cds_start0, cds_end0)

        motifs = []
        for left, right in zip(exons, exons[1:]):
            intron = genome.fetch(str(model["chrom"]), left[1], right[0]).upper()
            oriented = intron if model["strand"] == "+" else reverse_complement(intron)
            motifs.append(oriented[:2] + "-" + oriented[-2:] if len(oriented) >= 4 else "SHORT_INTRON")

        locus_id = f"R1M6_CAND{index:02d}"
        gene_start = min(start for start, _ in exons) + 1
        gene_end = max(end for _, end in exons)
        gff_lines.append(
            f'{model["chrom"]}\tR1M6_validation\tgene\t{gene_start}\t{gene_end}\t.\t{model["strand"]}\t.\tID={locus_id};status=provisional_not_in_catalogue\n'
        )
        gff_lines.append(
            f'{model["chrom"]}\tR1M6_validation\tmRNA\t{gene_start}\t{gene_end}\t.\t{model["strand"]}\t.\tID={locus_id}.t1;Parent={locus_id};source_transcript={tid}\n'
        )
        cumulative = 0
        for part_number, (start0, end0) in enumerate(projected, start=1):
            phase = (3 - cumulative % 3) % 3
            gff_lines.append(
                f'{model["chrom"]}\tR1M6_validation\tCDS\t{start0 + 1}\t{end0}\t.\t{model["strand"]}\t{phase}\tID={locus_id}.cds{part_number};Parent={locus_id}.t1\n'
            )
            cumulative += end0 - start0
        audit.append(
            {
                "provisional_locus_id": locus_id,
                "transcript_id": tid,
                "orf_id": oid,
                "chrom": model["chrom"],
                "model_start": gene_start,
                "model_end": gene_end,
                "strand": model["strand"],
                "exon_count": len(exons),
                "splice_motifs_transcript_orientation": ";".join(motifs) if motifs else "SINGLE_EXON",
                "all_splice_motifs_GT_AG": "YES" if motifs and all(value == "GT-AG" for value in motifs) else ("NA" if not motifs else "NO"),
                "cds_start_in_transcript_1based": cds_start0 + 1,
                "cds_end_in_transcript_1based": cds_end0,
                "downstream_stop_codon": stop,
                "protein_length_aa": len(peptide),
                "protein_sequence": peptide,
                "translated_sequence_exact_match": "YES",
                "genomic_cds_segments_transcript_order": ";".join(f"{start + 1}-{end}" for start, end in projected),
                "catalogue_status": "PROVISIONAL_NOT_ADDED_TO_CATALOGUE",
            }
        )

    write_path = out / "priority_candidate_genomic_model_audit.tsv"
    with write_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit)
    (out / "priority_candidate_provisional_models.gff3").write_text("".join(gff_lines), encoding="utf-8")
    print(f"Wrote {len(audit)} priority candidate model audits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
