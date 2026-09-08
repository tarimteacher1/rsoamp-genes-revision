#!/usr/bin/env python3
"""Apply the R. soongarica nsLTP non-phylogenetic gates to Tamarix candidates."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


SPECIFIC = {"IPR000528", "IPR033872", "IPR044741", "IPR039265"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preliminary", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--predgpi", required=True, type=Path)
    parser.add_argument("--tmbed", required=True, type=Path)
    parser.add_argument("--similarity", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--canonical-fasta", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


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
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    return records


def interval_end(interval: str) -> int | None:
    match = re.match(r"\d+-(\d+)", interval or "")
    return int(match.group(1)) if match else None


def eight_cm(sequence: str) -> tuple[str, str]:
    positions = [index + 1 for index, aa in enumerate(sequence) if aa == "C"]
    if len(positions) < 8:
        return "NO", ""
    best = None
    for offset in range(len(positions) - 7):
        subset = positions[offset : offset + 8]
        span = subset[-1] - subset[0]
        if best is None or span < best[0]:
            best = (span, subset)
    subset = best[1]
    spacing = "-".join(str(subset[index + 1] - subset[index] - 1) for index in range(7))
    return "YES", spacing


def main() -> None:
    args = parse_args()
    preliminary = read_tsv(args.preliminary, "protein_id")
    sequences = read_fasta(args.fasta)
    predgpi = read_tsv(args.predgpi, "protein_id")
    tmbed = read_tsv(args.tmbed, "protein_id")
    similarity = read_tsv(args.similarity, "protein_id")
    candidate_ids = [protein_id for protein_id, row in preliminary.items() if "nsLTP" in row["hmm_families"].split(";")]
    if set(candidate_ids) != set(sequences) or set(sequences) != set(predgpi) or set(predgpi) != set(tmbed) or set(tmbed) != set(similarity):
        raise ValueError("Tamarix nsLTP evidence inputs do not contain identical candidate IDs")

    rows = []
    for protein_id in sorted(candidate_ids):
        base = preliminary[protein_id]
        tm = tmbed[protein_id]
        sim = similarity[protein_id]
        gpi = predgpi[protein_id]
        deep_end = int(base["deepsig_cleavage_after"]) if base["deepsig_cleavage_after"] else None
        tm_end = interval_end(tm["signal_interval"]) if tm["signal_peptide"] == "YES" else None
        cleavage = deep_end or tm_end
        secreted = cleavage is not None
        mature = sequences[protein_id][cleavage:] if cleavage else sequences[protein_id]
        motif, spacing = eight_cm(mature)
        specific = sorted(SPECIFIC & set(base["interpro_accessions"].split(";")))
        margin = float(sim["positive_minus_strongest_alternative_bitscore"])
        positive = sim["reference_similarity_interpretation"] == "positive_reference_preferred"
        gpi_yes = gpi["predgpi_GPI_anchor"] == "YES"

        if secreted and specific and (motif == "YES" or (positive and margin >= 20)):
            classification = "A_CANONICAL_NSLTP"
            subtype = "LTPg_PREDICTED" if gpi_yes else "CANONICAL_CORE"
            confidence = "HIGH" if motif == "YES" and deep_end and tm_end else "MODERATE"
            reason = "Same-rule R1: secretion support, nsLTP-specific InterPro evidence, and compatible 8CM or strong curated homology"
        elif secreted and motif == "YES" and gpi_yes and positive and margin >= 50:
            classification = "A_CANONICAL_NSLTP"
            subtype = "LTPg_PREDICTED"
            confidence = "MODERATE"
            reason = "Same-rule R2: secretion, complete 8CM, PredGPI anchor, and curated-positive margin >=50 bits"
        else:
            classification = "C_EXCLUDED_OR_AMBIGUOUS_PROLAMIN_BOUNDARY"
            subtype = "NOT_COUNTED"
            confidence = "MODERATE"
            reason = "Did not meet the same canonical-core or LTPg evidence gate used for R. soongarica"

        rows.append(
            {
                **base,
                "tmbed_signal_peptide": tm["signal_peptide"],
                "tmbed_signal_interval": tm["signal_interval"],
                "tmbed_post_signal_tm_helix": tm["post_signal_tm_helix"],
                "combined_signal_peptide": "YES" if secreted else "NO",
                "effective_cleavage_after": cleavage or "",
                "recomputed_mature_cysteines": mature.count("C"),
                "recomputed_candidate_8CM": motif,
                "recomputed_8CM_spacing": spacing,
                "predgpi_GPI_anchor": gpi["predgpi_GPI_anchor"],
                "predgpi_omega_site": gpi["predgpi_omega_site"],
                "best_positive_reference": sim["best_positive_subject"],
                "best_positive_bitscore": sim["best_positive_bitscore"],
                "best_ambiguous_prolamin_reference": sim["best_ambiguous_prolamin_subject"],
                "best_ambiguous_prolamin_bitscore": sim["best_ambiguous_prolamin_bitscore"],
                "positive_minus_strongest_alternative_bitscore": sim["positive_minus_strongest_alternative_bitscore"],
                "reference_similarity_interpretation": sim["reference_similarity_interpretation"],
                "nsLTP_specific_IPR": ";".join(specific),
                "final_nsLTP_classification": classification,
                "final_nsLTP_subtype": subtype,
                "classification_confidence": confidence,
                "classification_reason": reason,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    canonical = [row for row in rows if row["final_nsLTP_classification"] == "A_CANONICAL_NSLTP"]
    with args.canonical_fasta.open("w", encoding="ascii") as handle:
        for row in canonical:
            handle.write(f">{row['protein_id']}\n{sequences[row['protein_id']]}\n")
    counts = Counter(row["final_nsLTP_classification"] for row in rows)
    count_path = args.output.parent / "tamarix_nsLTP_final_counts.tsv"
    with count_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["classification", "n"])
        writer.writerows(sorted(counts.items()))
    old_high = {protein_id for protein_id, row in preliminary.items() if row["final_family"] == "nsLTP" and row["classification_status"] == "HIGH_CONFIDENCE_SECRETED_FAMILY_MEMBER"}
    final_high = {row["protein_id"] for row in canonical}
    comparison = args.output.parent / "tamarix_nsLTP_preliminary_vs_final.tsv"
    with comparison.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_id", "preliminary_high_confidence", "final_canonical", "change"])
        for protein_id in sorted(old_high | final_high):
            writer.writerow([
                protein_id,
                "YES" if protein_id in old_high else "NO",
                "YES" if protein_id in final_high else "NO",
                "UNCHANGED" if (protein_id in old_high) == (protein_id in final_high) else "CHANGED",
            ])
    print(f"Tamarix same-rule nsLTP curation: candidates={len(rows)} canonical={len(canonical)} changed={len(old_high ^ final_high)}")


if __name__ == "__main__":
    main()
