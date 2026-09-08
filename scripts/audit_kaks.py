#!/usr/bin/env python3
"""Audit historical Ka/Ks pairs for catalogue membership and short-sequence reliability."""

from __future__ import annotations

import argparse
import csv
import math
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--catalogue", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


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
                    records[current] = "".join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks).upper()
    return records


def translate(cds: str) -> str:
    return "".join(CODON_TABLE.get(cds[i : i + 3], "X") for i in range(0, len(cds) - 2, 3)).rstrip("*")


def read_axt(path: Path) -> dict[str, tuple[str, str]]:
    blocks = {}
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    index = 0
    while index < len(lines):
        if not lines[index]:
            index += 1
            continue
        if index + 2 >= len(lines):
            raise ValueError("Incomplete AXT block")
        blocks[lines[index]] = (lines[index + 1], lines[index + 2])
        index += 3
    return blocks


def finite_float(value: str) -> float | None:
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_possibly_empty_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    catalogue = read_table(args.catalogue)
    included = {row["member_id"]: row for row in catalogue if row["catalogue_status"].startswith("INCLUDED")}
    historical = {row["member_id"]: row for row in read_table(args.project / "tables/T1_AMP_members.csv", ",")}
    pairs = read_table(args.project / "intermediate/dup_pairs.tsv")
    raw_results = {
        row["Sequence"]: row for row in read_table(args.project / "intermediate/kaks_out.txt")
    }
    axt = read_axt(args.project / "intermediate/kaks_input.axt")
    proteins = read_fasta(args.project / "data/protein.fa")
    cds = read_fasta(args.project / "data/cds.fa")

    rows = []
    for pair in pairs:
        member_a, member_b = pair["member_a"], pair["member_b"]
        pair_name = f"{member_a}-{member_b}"
        result = raw_results.get(pair_name, {})
        alignment = axt.get(pair_name)
        ka = finite_float(result.get("Ka", ""))
        ks = finite_float(result.get("Ks", ""))
        ratio = finite_float(result.get("Ka/Ks", ""))
        substitutions = int(float(result.get("Substitutions", 0) or 0))
        fisher_p = finite_float(result.get("P-Value(Fisher)", ""))
        alignment_nt = len(alignment[0]) if alignment else 0
        internal_stops = 0
        if alignment:
            internal_stops = sum(translate(sequence).count("*") for sequence in alignment)

        source_translation = []
        for member in (member_a, member_b):
            protein_id = historical[member]["protein_id"]
            if protein_id not in proteins or protein_id not in cds:
                source_translation.append(False)
            else:
                source_translation.append(translate(cds[protein_id]) == proteins[protein_id].rstrip("*"))

        flags = []
        if member_a not in included or member_b not in included:
            flags.append("PAIR_CONTAINS_MEMBER_EXCLUDED_FROM_FINAL_CATALOGUE")
        if not result:
            flags.append("NO_KAKS_RESULT")
        if alignment_nt < 150:
            flags.append("CODON_ALIGNMENT_LT_150_NT")
        if substitutions < 10:
            flags.append("FEWER_THAN_10_SUBSTITUTIONS")
        if ks is None or ks <= 0:
            flags.append("DS_NONPOSITIVE_OR_NONFINITE")
        elif ks > 2:
            flags.append("DS_GT_2_POTENTIAL_SATURATION")
        if ratio is None:
            flags.append("RATIO_NONFINITE")
        if internal_stops:
            flags.append("INTERNAL_STOP_IN_FILTERED_CODON_ALIGNMENT")
        if not all(source_translation):
            flags.append("SOURCE_CDS_TRANSLATION_MISMATCH")

        if "PAIR_CONTAINS_MEMBER_EXCLUDED_FROM_FINAL_CATALOGUE" in flags:
            reliability = "NOT_APPLICABLE_FINAL_CATALOGUE"
        elif any(flag in flags for flag in ("NO_KAKS_RESULT", "DS_NONPOSITIVE_OR_NONFINITE", "RATIO_NONFINITE", "INTERNAL_STOP_IN_FILTERED_CODON_ALIGNMENT", "SOURCE_CDS_TRANSLATION_MISMATCH")):
            reliability = "UNRELIABLE"
        elif "DS_GT_2_POTENTIAL_SATURATION" in flags or "FEWER_THAN_10_SUBSTITUTIONS" in flags or "CODON_ALIGNMENT_LT_150_NT" in flags:
            reliability = "UNRELIABLE_OR_HIGH_VARIANCE"
        else:
            reliability = "USABLE_WITH_SHORT_PEPTIDE_CAUTION"

        supported_below_one = (
            reliability == "USABLE_WITH_SHORT_PEPTIDE_CAUTION"
            and ratio is not None and ratio < 1
            and fisher_p is not None and fisher_p < 0.05
        )
        rows.append(
            {
                "pair": pair_name,
                "member_a": member_a,
                "member_b": member_b,
                "family": historical[member_a]["family"],
                "duplication_mode": pair["mode"],
                "method": result.get("Method", "YN" if result else ""),
                "Ka": result.get("Ka", ""),
                "Ks": result.get("Ks", ""),
                "Ka_Ks": result.get("Ka/Ks", ""),
                "fisher_p": result.get("P-Value(Fisher)", ""),
                "codon_alignment_length_nt": alignment_nt,
                "aligned_codons_after_complete_gap_deletion": alignment_nt // 3,
                "substitutions": substitutions,
                "source_cds_translate_exactly_to_protein": "YES" if all(source_translation) else "NO",
                "internal_stop_codons_after_complete_gap_deletion": internal_stops,
                "ratio_below_one": "YES" if ratio is not None and ratio < 1 else "NO" if ratio is not None else "",
                "fisher_supported_ratio_below_one": "YES" if supported_below_one else "NO",
                "quality_flags": ";".join(flags) if flags else "NONE",
                "reliability_class": reliability,
                "interpretation": "Consistent with purifying selection" if supported_below_one else "No robust pair-level purifying-selection inference",
            }
        )

    write_tsv(args.outdir / "kaks_quality_audit.tsv", rows)
    reliable = [row for row in rows if row["reliability_class"] == "USABLE_WITH_SHORT_PEPTIDE_CAUTION"]
    unreliable = [row for row in rows if row["reliability_class"] != "USABLE_WITH_SHORT_PEPTIDE_CAUTION"]
    fieldnames = list(rows[0])
    write_possibly_empty_tsv(args.outdir / "kaks_reliable_pairs.tsv", reliable, fieldnames)
    write_possibly_empty_tsv(args.outdir / "kaks_unreliable_pairs.tsv", unreliable, fieldnames)

    supported = sum(row["fisher_supported_ratio_below_one"] == "YES" for row in reliable)
    limitation = f"""# Ka/Ks reliability audit

Status: PASS

The historical analysis used protein-guided MAFFT alignment, complete deletion of codons containing alignment gaps, and the YN method in KaKs_Calculator. All {len(rows)} nominated duplicate pairs were retained in the audit, including pairs without a finite result or no longer applicable to the final catalogue.

After final-catalogue membership, source-CDS translation, alignment length, substitution count, finite dS, and dS-saturation checks, {len(reliable)} pairs were classified as usable with short-peptide caution. Of these, {supported} had Ka/Ks < 1 with Fisher P < 0.05 and can be described as *consistent with* purifying selection. Other ratios are descriptive only.

AMP coding sequences are short, so dN and dS estimates can have high sampling variance. dS > 2 was treated as potential saturation, and complete deletion of gap-containing codons can remove structurally variable regions. Ka/Ks therefore does not prove functional constraint or identify the biological source of selection.
"""
    (args.outdir / "kaks_limitation_text.md").write_text(limitation, encoding="utf-8")
    (args.outdir / "KAKS_AUDIT_COMPLETE.PASS").write_text(
        f"PASS\tnominated_pairs={len(rows)}\treliable_with_caution={len(reliable)}\tfisher_supported_below_one={supported}\n",
        encoding="utf-8",
    )
    print(limitation)


if __name__ == "__main__":
    main()
