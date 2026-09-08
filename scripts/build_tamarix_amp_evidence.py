#!/usr/bin/env python3
"""Build a traceable AMP-family evidence matrix for Tamarix austromongolica."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path


REVISION = Path("/path/to/rsoamp/revision_R1_20260902")
OUT = REVISION / "05_comparative_genomics"
PRED = OUT / "predictions"
INPUT = PRED / "tamarix_amp_hmm_candidates.faa"
INTERPRO = PRED / "tamarix_amp_interproscan.tsv"
DEEPSIG = PRED / "tamarix_amp_candidates.deepsig.gff3"
GFF = OUT / "inputs/Tamarix_austromongolica/tau.longest.gff3"

NSLTP_SPECIFIC_IPR = {"IPR000528", "IPR033872", "IPR044741", "IPR039265"}
FAMILY_SIGNATURES = {
    "Defensin": {"PF00304", "IPR008176"},
    "Snakin_GASA": {"PF02704", "IPR003854"},
    "Thionin": {"PF00321"},
    "Hevein_like": {"PF00187"},
    "Cyclotide": {"PF03784"},
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
            elif current is not None:
                records[current].append(line)
    return {key: "".join(value).replace("*", "") for key, value in records.items()}


def parse_deepsig(path: Path) -> dict[str, int]:
    cleavage = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[2] == "Signal peptide":
                cleavage[fields[0]] = int(fields[4])
    return cleavage


def parse_gff(path: Path) -> dict[str, tuple[str, int, int, str]]:
    coordinates = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "mRNA":
                continue
            attrs = dict(part.split("=", 1) for part in fields[8].split(";") if "=" in part)
            transcript = attrs.get("ID")
            if transcript:
                coordinates[transcript] = (fields[0], int(fields[3]), int(fields[4]), fields[6])
    return coordinates


def eight_cm(seq: str) -> tuple[str, str]:
    match = re.search(r"C([^C]{1,30})C([^C]{1,35})CC([^C]{1,35})C(.)C([^C]{1,35})C([^C]{1,35})C", seq)
    if not match:
        return "NO", ""
    cysteines = [index + 1 for index, aa in enumerate(match.group(0)) if aa == "C"]
    spacing = "-".join(str(b - a - 1) for a, b in zip(cysteines, cysteines[1:]))
    return "YES", spacing


def write_fasta(path: Path, rows: list[dict], sequences: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            protein_id = row["protein_id"]
            handle.write(f">TAU|{protein_id}|{row['final_family']}\n{sequences[protein_id]}\n")


def main() -> int:
    sequences = read_fasta(INPUT)
    cleavage = parse_deepsig(DEEPSIG)
    coordinates = parse_gff(GFF)

    hmm = defaultdict(list)
    with (PRED / "tamarix_amp_hmm_hits.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hmm[row["protein_id"]].append(row)

    interpro = defaultdict(list)
    with INTERPRO.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 13:
                interpro[fields[0]].append(fields)

    rows = []
    for protein_id, sequence in sequences.items():
        ipr_rows = interpro.get(protein_id, [])
        signatures = {row[4] for row in ipr_rows}
        ipr_accessions = {row[11] for row in ipr_rows if row[11] != "-"}
        descriptions = {row[12] for row in ipr_rows if row[12] != "-"}
        analyses = {row[3] for row in ipr_rows}
        sp_end = cleavage.get(protein_id, 0)
        mature = sequence[sp_end:]
        motif, spacing = eight_cm(mature)

        exact_families = [
            family
            for family, accessions in FAMILY_SIGNATURES.items()
            if accessions & (signatures | ipr_accessions)
        ]
        if NSLTP_SPECIFIC_IPR & ipr_accessions:
            exact_families.append("nsLTP")

        hmm_families = sorted({row["family_hmm"] for row in hmm.get(protein_id, [])})
        final_family = ""
        status = "REVIEW"
        reason = ""
        if len(set(exact_families)) > 1:
            final_family = "multi_domain_or_conflicting"
            status = "EXCLUDE_FROM_CORE_COMPARISON"
            reason = "Multiple family-specific signatures were detected."
        elif exact_families:
            final_family = exact_families[0]
            min_cys = 8 if final_family in {"Defensin", "nsLTP"} else 10 if final_family == "Snakin_GASA" else 6
            architecture_ok = mature.count("C") >= min_cys
            if sp_end and architecture_ok and (final_family != "nsLTP" or motif == "YES"):
                status = "HIGH_CONFIDENCE_SECRETED_FAMILY_MEMBER"
                reason = "Family-specific domain, N-terminal signal peptide and cysteine architecture agree."
            elif sp_end and architecture_ok:
                status = "FAMILY_SUPPORTED_ARCHITECTURE_REVIEW"
                reason = "Family-specific domain and secretion are supported, but the strict family motif screen failed."
            else:
                status = "FAMILY_SUPPORTED_ANNOTATION_REVIEW"
                reason = "Family-specific domain is present, but secretion or cysteine architecture is incomplete."
        elif "PF00234" in signatures or "nsLTP" in hmm_families:
            final_family = "nsLTP_prolamin_boundary"
            status = "AMBIGUOUS_NOT_COUNTED"
            reason = "Only broad PF00234/HMM evidence was detected without an nsLTP-specific InterPro term."
        elif "Hevein_like" in hmm_families:
            final_family = "chitin_binding_domain_protein"
            status = "EXCLUDE_FROM_CORE_COMPARISON"
            reason = "PF00187 is not sufficient to identify a hevein-like antimicrobial peptide."
        else:
            final_family = "unresolved_HMM_candidate"
            status = "AMBIGUOUS_NOT_COUNTED"
            reason = "No family-specific InterPro term was recovered."

        chrom, start, end, strand = coordinates.get(protein_id, ("", "", "", ""))
        rows.append(
            {
                "protein_id": protein_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "length_aa": len(sequence),
                "hmm_families": ";".join(hmm_families),
                "hmm_best_evalue": min((float(row["full_evalue"]) for row in hmm.get(protein_id, [])), default=""),
                "interpro_analyses": ";".join(sorted(analyses)),
                "signature_accessions": ";".join(sorted(signatures)),
                "interpro_accessions": ";".join(sorted(ipr_accessions)),
                "interpro_descriptions": ";".join(sorted(descriptions)),
                "deepsig_signal_peptide": "YES" if sp_end else "NO",
                "deepsig_cleavage_after": sp_end or "",
                "mature_cysteines": mature.count("C"),
                "candidate_8CM": motif,
                "candidate_8CM_spacing": spacing,
                "final_family": final_family,
                "classification_status": status,
                "classification_reason": reason,
            }
        )

    matrix = OUT / "tamarix_amp_evidence_preliminary.tsv"
    with matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    accepted = [row for row in rows if row["classification_status"] == "HIGH_CONFIDENCE_SECRETED_FAMILY_MEMBER"]
    counts = defaultdict(int)
    for row in rows:
        counts[(row["final_family"], row["classification_status"])] += 1
    with (OUT / "tamarix_amp_classification_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["final_family", "classification_status", "n_candidates"])
        for key in sorted(counts):
            writer.writerow([*key, counts[key]])

    for family in ("Defensin", "Snakin_GASA", "nsLTP"):
        write_fasta(
            PRED / f"tamarix_{family}_high_confidence.faa",
            [row for row in accepted if row["final_family"] == family],
            sequences,
        )
    print(f"Wrote {len(rows)} Tamarix candidates; {len(accepted)} pass the preliminary high-confidence gate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
