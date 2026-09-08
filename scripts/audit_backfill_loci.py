#!/usr/bin/env python3
"""Reconstruct the genomic evidence behind the six historical tblastn backfills."""

from __future__ import annotations

import csv
import subprocess
from collections import defaultdict
from pathlib import Path


ORIGINAL = Path("/path/to/rsoamp")
REVISION = ORIGINAL / "revision_R1_20260902"
OUT = REVISION / "03_annotation_rescue"
SAMTOOLS = Path("/path/to/user-home/tools/miniconda3/envs/as_splice/bin/samtools")
BACKFILL_IDS = {"RsDEF2", "RsDEF3", "RsDEF4", "RsDEF5", "RsDEF7", "RsSNA13"}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
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
            elif current:
                records[current].append(line)
    return {key: "".join(value).replace("*", "") for key, value in records.items()}


def fetch(chrom: str, start: int, end: int) -> str:
    result = subprocess.run(
        [str(SAMTOOLS), "faidx", str(ORIGINAL / "data/genome.fa"), f"{chrom}:{start}-{end}"],
        check=True,
        capture_output=True,
        text=True,
    )
    return "".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(">"))


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str, frame: int) -> str:
    seq = seq.upper()
    return "".join(CODON_TABLE.get(seq[i : i + 3], "X") for i in range(frame, len(seq) - 2, 3))


def parse_attributes(text: str) -> dict[str, str]:
    out = {}
    for part in text.split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
    return out


def load_gff() -> dict[str, list[dict]]:
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    with (ORIGINAL / "data/genome.gff").open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            by_chrom[fields[0]].append(
                {
                    "type": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "attributes": parse_attributes(fields[8]),
                }
            )
    return by_chrom


def nearest_gene(features: list[dict], start: int, end: int) -> tuple[str, int | str, str]:
    genes = [feature for feature in features if feature["type"] == "gene"]
    if not genes:
        return "", "", ""
    def distance(feature: dict) -> int:
        if feature["end"] < start:
            return start - feature["end"]
        if feature["start"] > end:
            return feature["start"] - end
        return 0
    gene = min(genes, key=distance)
    return gene["attributes"].get("ID", ""), distance(gene), f"{gene['start']}-{gene['end']}({gene['strand']})"


def locate_peptide(chrom: str, locus_start: int, locus_end: int, strand: str, peptide: str) -> dict:
    window_start = max(1, locus_start - 1500)
    window_end = locus_end + 1500
    genomic = fetch(chrom, window_start, window_end).upper()
    oriented = genomic if strand == "+" else revcomp(genomic)
    matches = []
    for frame in range(3):
        aa = translate(oriented, frame)
        offset = aa.find(peptide)
        while offset >= 0:
            matches.append((frame, offset, aa))
            offset = aa.find(peptide, offset + 1)
    if len(matches) != 1:
        return {
            "exact_genomic_encoding": "NO" if not matches else "MULTIPLE",
            "exact_match_count": len(matches),
            "coding_start": "",
            "coding_end": "",
            "frame": "",
            "upstream_inframe_methionine": "",
            "amino_acids_from_nearest_methionine": "",
            "downstream_stop_before_next_methionine": "",
            "next_amino_acid": "",
        }
    frame, offset, aa = matches[0]
    nt_offset = frame + offset * 3
    if strand == "+":
        coding_start = window_start + nt_offset
        coding_end = coding_start + len(peptide) * 3 - 1
    else:
        coding_end = window_end - nt_offset
        coding_start = coding_end - len(peptide) * 3 + 1

    left_stop = aa.rfind("*", 0, offset)
    orf_prefix = aa[left_stop + 1 : offset]
    methionine = orf_prefix.rfind("M")
    from_methionine = len(orf_prefix) - methionine if methionine >= 0 else ""
    after = offset + len(peptide)
    right_stop = aa.find("*", after)
    next_m = aa.find("M", after, right_stop if right_stop >= 0 else len(aa))
    return {
        "exact_genomic_encoding": "YES",
        "exact_match_count": 1,
        "coding_start": coding_start,
        "coding_end": coding_end,
        "frame": frame,
        "upstream_inframe_methionine": "YES" if methionine >= 0 else "NO",
        "amino_acids_from_nearest_methionine": from_methionine,
        "downstream_stop_before_next_methionine": "YES" if right_stop >= 0 and (next_m < 0 or right_stop < next_m) else "NO",
        "next_amino_acid": aa[after] if after < len(aa) else "",
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    members = [
        row for row in csv.DictReader((ORIGINAL / "tables/T1_AMP_members.csv").open(encoding="utf-8-sig"))
        if row["member_id"] in BACKFILL_IDS
    ]
    peptides = read_fasta(ORIGINAL / "intermediate/members_final.faa")
    gff = load_gff()
    output = []
    for member in members:
        peptide = peptides[member["member_id"]]
        chrom = member["chrom"]
        start, end = int(member["start"]), int(member["end"])
        located = locate_peptide(chrom, start, end, member["strand"], peptide)
        feature_start = int(located["coding_start"] or start)
        feature_end = int(located["coding_end"] or end)
        overlapping = [
            feature for feature in gff.get(chrom, [])
            if feature["start"] <= feature_end and feature["end"] >= feature_start
        ]
        gene_id, gene_distance, gene_coords = nearest_gene(gff.get(chrom, []), feature_start, feature_end)
        local_seq = fetch(chrom, max(1, feature_start - 2000), feature_end + 2000).upper()
        output.append(
            {
                "member_id": member["member_id"],
                "family": member["family"],
                "historical_locus": f"{chrom}:{start}-{end}({member['strand']})",
                "historical_locus_nt": end - start + 1,
                "hardcoded_peptide_length_aa": len(peptide),
                **located,
                "exact_coding_interval_within_historical_locus": (
                    "YES" if located["coding_start"] and start <= int(located["coding_start"]) <= int(located["coding_end"]) <= end else "NO"
                ),
                "overlapping_gff_feature_types": ";".join(sorted({feature["type"] for feature in overlapping})),
                "overlapping_gff_ids": ";".join(sorted({feature["attributes"].get("ID", "") for feature in overlapping if feature["attributes"].get("ID")})),
                "nearest_gene_id": gene_id,
                "nearest_gene_distance_bp": gene_distance,
                "nearest_gene_coordinates": gene_coords,
                "N_bases_within_2kb_flanks": local_seq.count("N"),
                "rna_short_read_support": "PENDING_VERIFIED_FASTQ_ALIGNMENT",
                "long_read_transcript_support": "NOT_RUN",
                "final_gene_model_status": "PENDING_ORF_SPLICE_AND_READ_SUPPORT_REVIEW",
            }
        )

    matrix = OUT / "backfill_locus_audit_preliminary.tsv"
    with matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)

    summary = OUT / "backfill_locus_audit_preliminary.md"
    exact = sum(row["exact_genomic_encoding"] == "YES" for row in output)
    with_start = sum(row["upstream_inframe_methionine"] == "YES" for row in output)
    summary.write_text(
        "# Six historical tblastn backfills - preliminary audit\n\n"
        "Status: **PARTIAL; no locus is accepted as a final gene model at this stage**\n\n"
        f"The six hardcoded peptide sequences were tested against the reported genomic loci and strands. "
        f"Exact genomic encoding was found for {exact}/6 peptides, and {with_start}/6 had an upstream in-frame methionine "
        "within the local stop-delimited interval. The matrix also records GFF overlap, nearest genes and assembly-gap screening.\n\n"
        "The old workflow used tblastn-HSP spans and absence of overlap with a GFF gene feature, then assumed a single exon. "
        "That procedure does not establish an intact gene model. Final calls require explicit ORF boundaries, splice/read evidence, "
        "domain and secretion evidence, and a genome-wide rescue search under the frozen family rules.\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(output)} rows to {matrix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
