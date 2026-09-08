#!/usr/bin/env python3
"""Audit AMP promoters and test motif-presence enrichment against genomic promoters."""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


PROMOTER_LENGTH = 2000
ELEMENTS = {
    "ABRE": ("ACGTG", "ABA-responsive"),
    "G-box": ("CACGTG", "light/ABA"),
    "DRE_CRT": ("RCCGAC", "dehydration/cold"),
    "MYB": ("CNGTTR", "drought-associated MYB consensus"),
    "MYB_MBS": ("CAACTG", "drought-associated MBS consensus"),
    "MYC": ("CANNTG", "MYC consensus"),
    "W-box": ("TTGACY", "WRKY-associated W-box consensus"),
    "TC-rich": ("GTTTTCTTAC", "stress-associated consensus"),
    "LTR": ("CCGAAA", "low-temperature-associated consensus"),
    "GT-1": ("GAAAAA", "GT-1-associated consensus"),
    "as-1_TGACG": ("TGACG", "as-1-like consensus"),
    "CGTCA-motif": ("CGTCA", "methyl-jasmonate-associated consensus"),
    "ARE": ("AAACCA", "anaerobic-response-associated consensus"),
    "P-box_GARE": ("CCTTTTG", "gibberellin-associated consensus"),
    "ABRE-like": ("BACGTG", "ABRE-like consensus"),
    "STRE": ("AGGGG", "general-stress-associated consensus"),
}
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
IUPAC_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attributes[key] = value
    return attributes


def load_genes(path: Path) -> dict[str, dict[str, object]]:
    genes: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            # The frozen AMP catalogue uses the GFF gene Name (gene_*), whereas
            # ID stores the GeMoMa model identifier (gmmChr*).  Use Name as the
            # unique gene key so targets and the genome-wide background share
            # the same identifier namespace.
            gene_id = attrs.get("Name") or attrs.get("gene_id") or attrs.get("ID")
            if not gene_id:
                continue
            genes[gene_id] = {
                "gene_id": gene_id,
                "chrom": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            }
    if not genes:
        raise ValueError("No gene features were parsed from the GFF")
    return genes


class FastaIndex:
    def __init__(self, fasta: Path):
        self.fasta = fasta
        fai = Path(str(fasta) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"Missing FASTA index: {fai}")
        self.entries: dict[str, tuple[int, int, int, int]] = {}
        with fai.open(encoding="utf-8") as handle:
            for line in handle:
                name, length, offset, line_bases, line_width, *_ = line.rstrip("\n").split("\t")
                self.entries[name] = tuple(map(int, (length, offset, line_bases, line_width)))

    def length(self, chrom: str) -> int:
        return self.entries[chrom][0]

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        length, offset, line_bases, line_width = self.entries[chrom]
        start0 = max(0, start0)
        end0 = min(length, end0)
        if end0 <= start0:
            return ""
        byte_start = offset + (start0 // line_bases) * line_width + (start0 % line_bases)
        last = end0 - 1
        byte_end = offset + (last // line_bases) * line_width + (last % line_bases) + 1
        with self.fasta.open("rb") as handle:
            handle.seek(byte_start)
            raw = handle.read(byte_end - byte_start)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()[: end0 - start0]


def reverse_complement(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1]


def reverse_complement_iupac(consensus: str) -> str:
    return consensus.translate(IUPAC_COMPLEMENT)[::-1]


def motif_regex(consensus: str) -> re.Pattern[str]:
    return re.compile("(?=(" + "".join(IUPAC[base] for base in consensus) + "))")


def motif_count_both_strands(seq: str, consensus: str) -> int:
    sites: set[tuple[int, int]] = set()
    for motif in {consensus, reverse_complement_iupac(consensus)}:
        regex = motif_regex(motif)
        motif_length = len(motif)
        sites.update((match.start(), match.start() + motif_length) for match in regex.finditer(seq))
    return len(sites)


def promoter_interval(gene: dict[str, object], chromosome_length: int) -> tuple[int, int]:
    if gene["strand"] == "+":
        return max(0, int(gene["start"]) - 1 - PROMOTER_LENGTH), int(gene["start"]) - 1
    return int(gene["end"]), min(chromosome_length, int(gene["end"]) + PROMOTER_LENGTH)


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_greater(a: int, b: int, c: int, d: int) -> float:
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    maximum = min(row1, col1)
    terms = []
    for observed in range(a, maximum + 1):
        log_p = log_choose(col1, observed) + log_choose(total - col1, row1 - observed) - log_choose(total, row1)
        terms.append(math.exp(log_p))
    return min(1.0, sum(terms))


def bh_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 1.0
    n = len(p_values)
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * n / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    genes = load_genes(args.gff)
    genome = FastaIndex(args.genome)

    with args.catalogue.open(encoding="utf-8", newline="") as handle:
        catalogue = list(csv.DictReader(handle, delimiter="\t"))
    target_rows = [row for row in catalogue if row.get("gene_id") and row.get("catalogue_status", "INCLUDED").startswith("INCLUDED")]
    target_gene_ids = {row["gene_id"] for row in target_rows}
    missing = sorted(target_gene_ids - genes.keys())
    if missing:
        raise ValueError(f"Catalogue gene IDs absent from GFF: {', '.join(missing[:10])}")

    intervals_by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    starts_by_chrom: dict[str, list[int]] = {}
    for gene in genes.values():
        intervals_by_chrom[str(gene["chrom"])].append((int(gene["start"]), int(gene["end"]), str(gene["gene_id"])))
    for chrom in intervals_by_chrom:
        intervals_by_chrom[chrom].sort()
        starts_by_chrom[chrom] = [row[0] for row in intervals_by_chrom[chrom]]

    gene_motifs: dict[str, dict[str, int]] = {}
    audit_by_gene: dict[str, dict[str, object]] = {}
    for gene_id, gene in genes.items():
        chrom = str(gene["chrom"])
        if chrom not in genome.entries:
            continue
        start0, end0 = promoter_interval(gene, genome.length(chrom))
        sequence = genome.fetch(chrom, start0, end0)
        if gene["strand"] == "-":
            sequence = reverse_complement(sequence)
        counts = {name: motif_count_both_strands(sequence, consensus) for name, (consensus, _) in ELEMENTS.items()}
        gene_motifs[gene_id] = counts

        overlaps = []
        possible = intervals_by_chrom[chrom]
        cutoff = bisect.bisect_right(starts_by_chrom[chrom], end0)
        for other_start, other_end, other_id in possible[:cutoff]:
            if other_id != gene_id and other_end >= start0 + 1 and other_start <= end0:
                overlaps.append(other_id)
        audit_by_gene[gene_id] = {
            "gene_id": gene_id,
            "chrom": chrom,
            "gene_start": gene["start"],
            "gene_end": gene["end"],
            "strand": gene["strand"],
            "promoter_start_1based": start0 + 1 if end0 > start0 else "",
            "promoter_end_1based": end0,
            "promoter_length_bp": len(sequence),
            "scaffold_boundary_truncated": "YES" if len(sequence) < PROMOTER_LENGTH else "NO",
            "overlapping_other_gene_count": len(set(overlaps)),
            "overlapping_other_gene_ids": ";".join(sorted(set(overlaps))),
        }

    member_by_gene = {row["gene_id"]: row for row in target_rows}
    audit_rows: list[dict[str, object]] = []
    full_rows: list[dict[str, object]] = []
    for gene_id in sorted(target_gene_ids):
        member = member_by_gene[gene_id]
        audit_rows.append({"member_id": member.get("member_id", ""), "family": member.get("family", ""), **audit_by_gene[gene_id]})
        for element, (consensus, category) in ELEMENTS.items():
            count = gene_motifs[gene_id][element]
            full_rows.append(
                {
                    "member_id": member.get("member_id", ""),
                    "gene_id": gene_id,
                    "family": member.get("family", ""),
                    "element": element,
                    "iupac_consensus": consensus,
                    "annotation_category": category,
                    "nonredundant_genomic_site_count": count,
                    "present": "YES" if count > 0 else "NO",
                    "interpretation": "Sequence-level motif occurrence; no functional activation is inferred",
                }
            )

    target_full = [gene_id for gene_id in target_gene_ids if audit_by_gene[gene_id]["promoter_length_bp"] == PROMOTER_LENGTH]
    background_full = [
        gene_id for gene_id in gene_motifs
        if gene_id not in target_gene_ids and audit_by_gene[gene_id]["promoter_length_bp"] == PROMOTER_LENGTH
    ]
    enrichment_rows: list[dict[str, object]] = []
    p_values: list[float] = []
    for element, (consensus, category) in ELEMENTS.items():
        target_present = sum(gene_motifs[gene_id][element] > 0 for gene_id in target_full)
        background_present = sum(gene_motifs[gene_id][element] > 0 for gene_id in background_full)
        a = target_present
        b = len(target_full) - target_present
        c = background_present
        d = len(background_full) - background_present
        p_value = fisher_greater(a, b, c, d)
        p_values.append(p_value)
        enrichment_rows.append(
            {
                "element": element,
                "iupac_consensus": consensus,
                "annotation_category": category,
                "target_full_length_promoters": len(target_full),
                "target_with_motif": a,
                "target_fraction": f"{a / len(target_full):.6f}" if target_full else "",
                "non_target_background_full_length_promoters": len(background_full),
                "background_with_motif": c,
                "background_fraction": f"{c / len(background_full):.6f}" if background_full else "",
                "fisher_greater_p": f"{p_value:.8g}",
                "bh_adjusted_p": "",
                "enriched_at_bh_fdr_0.05": "",
            }
        )
    adjusted = bh_adjust(p_values)
    for row, value in zip(enrichment_rows, adjusted):
        row["bh_adjusted_p"] = f"{value:.8g}"
        row["enriched_at_bh_fdr_0.05"] = "YES" if value < 0.05 else "NO"

    write_tsv(args.outdir / "promoter_extraction_audit.tsv", audit_rows)
    write_tsv(args.outdir / "cis_element_full_results.tsv", full_rows)
    write_tsv(args.outdir / "cis_element_enrichment.tsv", enrichment_rows)

    significant = [row["element"] for row in enrichment_rows if row["enriched_at_bh_fdr_0.05"] == "YES"]
    interpretation = [
        "# Promoter motif audit",
        "",
        "Status: PASS",
        "",
        f"Promoters were defined as the {PROMOTER_LENGTH}-bp interval immediately upstream of each annotated gene in transcript orientation.",
        "Overlapping IUPAC matches were allowed, while forward- and reverse-complement matches at the same genomic coordinates were counted once.",
        "Scaffold-truncated promoters were retained in the occurrence table but excluded from enrichment testing.",
        f"The enrichment background comprised {len(background_full):,} complete promoters from non-target annotated genes; {len(target_full)} complete target promoters were tested.",
        "One-sided Fisher exact tests evaluated greater motif presence in target promoters, with Benjamini-Hochberg correction across the 16 tested motifs.",
        "",
        "Significant enriched motifs at BH FDR < 0.05: " + (", ".join(significant) if significant else "none"),
        "",
        "Motif occurrence or enrichment is a sequence-level association and does not establish transcription-factor binding, treatment-specific activation, or a regulatory mechanism.",
    ]
    (args.outdir / "promoter_interpretation.md").write_text("\n".join(interpretation) + "\n", encoding="utf-8")
    (args.outdir / "PROMOTER_AUDIT_COMPLETE.PASS").write_text(
        f"PASS\ttargets={len(target_gene_ids)}\ttarget_full={len(target_full)}\tbackground_full={len(background_full)}\n",
        encoding="utf-8",
    )
    print((args.outdir / "promoter_interpretation.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
