#!/usr/bin/env python3
"""Build a non-final nsLTP candidate evidence matrix without changing the old catalogue."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


ORIGINAL = Path("/path/to/rsoamp")
REVISION = ORIGINAL / "revision_R1_20260902"
OUT = REVISION / "01_nslTP_curation"
INPUT_FASTA = OUT / "inputs/nsLTP_candidate_pool.faa"
DEEPSIG_GFF = OUT / "predictions/nsLTP_candidates.deepsig.gff3"

NSLTP_SPECIFIC_IPR = {"IPR000528", "IPR033872", "IPR044741", "IPR039265"}
BROAD_PROLAMIN_TERMS = (
    "bifunctional inhibitor/plant lipid transfer protein/seed storage helical domain",
    "bifunctional inhibitor/lipid-transfer protein/seed storage 2s albumin superfamily protein",
    "alpha-amylase inhibitors (aai), lipid transfer (lt) and seed storage (ss) protein",
    "aai_ltss",
)

# Generic prolamin-superfamily annotations are not evidence for a competing
# function. Only specific, non-nsLTP calls belong in this field.
SPECIFIC_COMPETING_PATTERNS = (
    re.compile(r"(^|[^a-z0-9])2s (albumin|seed storage protein)([^a-z0-9]|$)"),
    re.compile(r"(^|[^a-z0-9])(alpha-)?amylase inhibitor([^a-z0-9]|$)"),
    re.compile(r"(^|[^a-z0-9])trypsin inhibitor([^a-z0-9]|$)"),
    re.compile(r"(^|[^a-z0-9])protease inhibitor([^a-z0-9]|$)"),
)


def is_broad_prolamin_annotation(text: str) -> bool:
    lowered = text.lower()
    return any(term in lowered for term in BROAD_PROLAMIN_TERMS)


def is_specific_competing_annotation(text: str) -> bool:
    lowered = text.lower()
    if is_broad_prolamin_annotation(lowered):
        return False
    return any(pattern.search(lowered) for pattern in SPECIFIC_COMPETING_PATTERNS)


def read_table(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_optional_tsv(path: Path, key: str) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {row[key]: row for row in read_table(path, "\t")}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line)
    return {key: "".join(value).replace("*", "") for key, value in records.items()}


def write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for name in sorted(records):
            handle.write(f">{name}\n")
            seq = records[name]
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")


def prepare() -> None:
    union = read_table(ORIGINAL / "intermediate/candidates_union.tsv", "\t")
    sequences = read_fasta(ORIGINAL / "intermediate/candidates.faa")
    candidate_ids = sorted({row["protein_id"] for row in union if row["family"] == "nsLTP"})
    missing = [protein_id for protein_id in candidate_ids if protein_id not in sequences]
    if missing:
        raise SystemExit(f"Sequences missing for {len(missing)} candidates: {missing}")
    write_fasta(INPUT_FASTA, {protein_id: sequences[protein_id] for protein_id in candidate_ids})
    print(f"Prepared {len(candidate_ids)} nsLTP/prolamin-superfamily candidates: {INPUT_FASTA}")


def parse_deepsig(path: Path) -> dict[str, int]:
    cleavage = {}
    if not path.exists():
        return cleavage
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[2] == "Signal peptide":
                cleavage[fields[0]] = int(fields[4])
    return cleavage


def eight_cm(seq: str) -> tuple[str, str, int, int]:
    # The screening expression records a candidate 8CM; final calls also require
    # curated alignment and domain/phylogeny evidence.
    pattern = re.compile(r"C[^C]{1,30}C[^C]{1,35}CC[^C]{1,35}C.C[^C]{1,35}C[^C]{1,35}C")
    match = pattern.search(seq)
    if not match:
        return "NO", "", 0, 0
    positions = [index + 1 for index, aa in enumerate(match.group()) if aa == "C"]
    spacing = ",".join(str(b - a - 1) for a, b in zip(positions, positions[1:]))
    return "YES", spacing, match.start() + 1, match.end()


def hydrophobic_tail_screen(seq: str, window: int = 19, threshold: float = 1.6) -> tuple[str, str, str]:
    kd = {
        "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
        "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
        "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
        "K": -3.9, "R": -4.5,
    }
    hits = []
    for start in range(35, max(35, len(seq) - window + 1)):
        segment = seq[start : start + window]
        if len(segment) != window or any(aa not in kd for aa in segment):
            continue
        score = sum(kd[aa] for aa in segment) / window
        if score >= threshold:
            hits.append((score, start + 1, start + window, segment))
    if not hits:
        return "NO", "", ""
    score, start, end, segment = max(hits)
    c_terminal = "YES" if end >= len(seq) - 10 else "NO"
    return c_terminal, f"{start}-{end}", f"{score:.3f}:{segment}"


def interval_end(interval: str) -> int | None:
    if not interval or "-" not in interval:
        return None
    try:
        return int(interval.rsplit("-", 1)[1])
    except ValueError:
        return None


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    union_rows = read_table(ORIGINAL / "intermediate/candidates_union.tsv", "\t")
    union = {row["protein_id"]: row for row in union_rows if row["family"] == "nsLTP"}
    sequences = read_fasta(INPUT_FASTA)
    t1 = read_table(ORIGINAL / "tables/T1_AMP_members.csv")
    current = {row["protein_id"]: row for row in t1 if row["family"] == "nsLTP"}
    old_t10 = {row["protein_id"]: row for row in read_table(ORIGINAL / "tables/T10_nsLTP_curation.csv")}
    rejected = {row["protein_id"]: row for row in read_table(ORIGINAL / "intermediate/rejected_candidates.csv") if row["prior_family"] == "nsLTP"}
    cleavage = parse_deepsig(DEEPSIG_GFF)
    reference_similarity = read_optional_tsv(OUT / "nsLTP_labeled_reference_similarity.tsv", "protein_id")
    predgpi = read_optional_tsv(OUT / "nsLTP_candidates_predgpi.tsv", "protein_id")
    tmbed = read_optional_tsv(OUT / "nsLTP_candidates_tmbed.tsv", "protein_id")

    ips = defaultdict(list)
    with (ORIGINAL / "intermediate/ips/candidates_ips.tsv").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 13 and fields[0] in union:
                ips[fields[0]].append(
                    {
                        "app": fields[3],
                        "signature": fields[4],
                        "description": fields[5],
                        "start": fields[6],
                        "end": fields[7],
                        "evalue": fields[8],
                        "ipr": "" if fields[11] == "-" else fields[11],
                        "ipr_name": "" if fields[12] == "-" else fields[12],
                    }
                )

    output = []
    for protein_id in sorted(union):
        seq = sequences[protein_id]
        sp_end = cleavage.get(protein_id)
        tm = tmbed.get(protein_id, {})
        tm_sp_end = interval_end(tm.get("signal_interval", "")) if tm.get("signal_peptide") == "YES" else None
        effective_sp_end = sp_end or tm_sp_end
        mature = seq[effective_sp_end:] if effective_sp_end else seq
        motif, spacing, motif_start, motif_end = eight_cm(mature)
        ctail, hydro_pos, hydro_detail = hydrophobic_tail_screen(seq)
        annotations = ips.get(protein_id, [])
        ipr_ids = sorted({row["ipr"] for row in annotations if row["ipr"]})
        signatures = sorted({f'{row["app"]}:{row["signature"]}' for row in annotations})
        descriptions = sorted({row["ipr_name"] or row["description"] for row in annotations})
        specific = sorted(set(ipr_ids) & NSLTP_SPECIFIC_IPR)
        broad_prolamin = sorted({text for text in descriptions if is_broad_prolamin_annotation(text)})
        competing = sorted({text for text in descriptions if is_specific_competing_annotation(text)})
        similarity = reference_similarity.get(protein_id, {})
        similarity_call = similarity.get("reference_similarity_interpretation", "NOT_RUN")
        gpi = predgpi.get(protein_id, {})
        gpi_call = gpi.get("predgpi_GPI_anchor", "NOT_RUN")
        post_signal_tm = tm.get("post_signal_tm_helix", "NOT_RUN")

        if sp_end and tm_sp_end:
            signal_concordance = "CONCORDANT_BOTH_POSITIVE"
        elif sp_end:
            signal_concordance = "DEEPSIG_ONLY"
        elif tm_sp_end:
            signal_concordance = "TMBED_ONLY"
        else:
            signal_concordance = "CONCORDANT_BOTH_NEGATIVE"

        if effective_sp_end and motif == "YES" and gpi_call == "YES" and post_signal_tm in {"YES", "NO"}:
            tier = "predicted_LTPg_supported_pending_phylogeny"
        elif effective_sp_end and motif == "YES" and post_signal_tm == "YES" and gpi_call != "YES":
            tier = "post_signal_transmembrane_candidate_pending"
        elif effective_sp_end and motif == "YES" and similarity_call == "broad_prolamin_umbrella_reference_preferred" and not specific:
            tier = "ambiguous_nsLTP_prolamin_boundary"
        elif effective_sp_end and motif == "YES" and similarity_call == "explicit_2S_seed_storage_reference_preferred" and not specific:
            tier = "likely_non_nsLTP_2S_pending_phylogeny"
        elif effective_sp_end and motif == "YES" and specific and not competing:
            tier = "strong_nsLTP_core_evidence"
        elif effective_sp_end and motif == "YES" and similarity_call == "positive_reference_preferred" and not competing:
            tier = "candidate_nsLTP_with_curated_homology"
        elif effective_sp_end and motif == "YES" and not competing:
            tier = "possible_nsLTP_core"
        elif competing:
            tier = "ambiguous_or_competing_prolamin"
        elif not effective_sp_end:
            tier = "partial_or_nonsecretory_pending"
        else:
            tier = "ambiguous_pending"

        output.append(
            {
                "member_id": current.get(protein_id, {}).get("member_id", ""),
                "protein_id": protein_id,
                "old_catalogue_status": old_t10.get(protein_id, {}).get("subclass", "excluded_or_not_carried_forward"),
                "old_rejection_reason": rejected.get(protein_id, {}).get("reason", ""),
                "length_aa": len(seq),
                "deepsig_signal_peptide": "YES" if sp_end else "NO",
                "deepsig_cleavage_after": sp_end or "",
                "combined_signal_peptide": "YES" if effective_sp_end else "NO",
                "effective_cleavage_after": effective_sp_end or "",
                "signal_prediction_concordance": signal_concordance,
                "estimated_mature_length": len(mature),
                "total_cysteines": seq.count("C"),
                "mature_cysteines": mature.count("C"),
                "candidate_8CM": motif,
                "candidate_8CM_spacing": spacing,
                "candidate_8CM_mature_coordinates": f"{motif_start}-{motif_end}" if motif == "YES" else "",
                "cterminal_hydrophobic_screen": ctail,
                "hydrophobic_window_coordinates": hydro_pos,
                "hydrophobic_window_detail": hydro_detail,
                "predgpi_GPI_anchor": gpi_call,
                "predgpi_omega_site": gpi.get("predgpi_omega_site", ""),
                "predgpi_score": gpi.get("predgpi_score", ""),
                "tmbed_signal_peptide": tm.get("signal_peptide", "NOT_RUN"),
                "tmbed_signal_interval": tm.get("signal_interval", ""),
                "tmbed_tm_alpha_helix_count": tm.get("tm_alpha_helix_count", ""),
                "tmbed_tm_alpha_intervals": tm.get("tm_alpha_intervals", ""),
                "tmbed_post_signal_tm_helix": post_signal_tm,
                "tmbed_post_signal_tm_intervals": tm.get("post_signal_tm_intervals", ""),
                "PF00234_evalue": union[protein_id]["hmm_evalue"],
                "PF00234_score": union[protein_id]["hmm_score"],
                "blast_seed": union[protein_id]["blast_seed"],
                "blast_family": union[protein_id]["blast_fam"],
                "blast_pident": union[protein_id]["blast_pident"],
                "blast_evalue": union[protein_id]["blast_evalue"],
                "nsLTP_specific_IPR": ";".join(specific),
                "all_IPR": ";".join(ipr_ids),
                "all_signatures": ";".join(signatures),
                "broad_prolamin_umbrella": ";".join(broad_prolamin),
                "competing_annotations": ";".join(competing),
                "best_positive_reference": similarity.get("best_positive_subject", ""),
                "best_positive_bitscore": similarity.get("best_positive_bitscore", ""),
                "best_explicit_2S_reference": similarity.get("best_explicit_2S_subject", ""),
                "best_explicit_2S_bitscore": similarity.get("best_explicit_2S_bitscore", ""),
                "best_ambiguous_prolamin_reference": similarity.get("best_ambiguous_prolamin_subject", ""),
                "best_ambiguous_prolamin_bitscore": similarity.get("best_ambiguous_prolamin_bitscore", ""),
                "positive_minus_strongest_alternative_bitscore": similarity.get("positive_minus_strongest_alternative_bitscore", ""),
                "reference_similarity_interpretation": similarity_call,
                "annotation_descriptions": ";".join(descriptions),
                "preliminary_evidence_tier": tier,
                "final_classification": "PENDING_REFERENCE_ALIGNMENT_TM_GPI_AND_PHYLOGENY",
                "classification_confidence": "PENDING",
            }
        )

    columns = list(output[0])
    matrix = OUT / "nsLTP_evidence_matrix_preliminary.tsv"
    with matrix.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)

    rsltp20 = next(row for row in output if row["member_id"] == "RsLTP20")
    report = f"""# RsLTP20 case report - preliminary evidence stage

Status: **OPEN; no final reclassification yet**

## Historical inconsistency

The old InterPro filtering file rejected `{rsltp20['protein_id']}` as having no AMP-specific domain, while the later script classified it as canonical whenever either an nsLTP-specific InterPro term **or any nsLTP-seed BLAST hit** was present. RsLTP20 was therefore restored solely through the BLAST branch of an overly permissive OR rule.

## Evidence currently observed

- Precursor length: {rsltp20['length_aa']} aa.
- DeepSig N-terminal signal peptide: {rsltp20['deepsig_signal_peptide']}; cleavage after residue {rsltp20['deepsig_cleavage_after']}.
- Candidate mature-region 8CM: {rsltp20['candidate_8CM']}; spacing `{rsltp20['candidate_8CM_spacing']}`.
- nsLTP-specific InterPro terms: `{rsltp20['nsLTP_specific_IPR'] or 'none'}`.
- Broad/other InterPro terms: `{rsltp20['all_IPR']}`.
- Best historical nsLTP-seed BLAST: `{rsltp20['blast_seed']}`, identity {rsltp20['blast_pident']}%, E-value {rsltp20['blast_evalue']}.
- C-terminal hydrophobic screening flag: {rsltp20['cterminal_hydrophobic_screen']}.
- PredGPI GPI-anchor prediction: {rsltp20['predgpi_GPI_anchor']}.
- TMbed post-signal transmembrane helix: {rsltp20['tmbed_post_signal_tm_helix']} ({rsltp20['tmbed_post_signal_tm_intervals'] or 'none'}).
- Best reviewed Arabidopsis nsLTP reference: `{rsltp20['best_positive_reference']}` (bitscore {rsltp20['best_positive_bitscore']}).
- Best broad prolamin-umbrella reference: `{rsltp20['best_ambiguous_prolamin_reference']}` (bitscore {rsltp20['best_ambiguous_prolamin_bitscore']}).
- Positive-minus-strongest-alternative bitscore: {rsltp20['positive_minus_strongest_alternative_bitscore']}.
- Labeled-reference interpretation: `{rsltp20['reference_similarity_interpretation']}`.

## Decision still required

The signal peptide and cysteine framework are compatible with an nsLTP hypothesis, but the old best hit was an ambiguous protease-inhibitor/seed-storage/LTP-superfamily sequence rather than a curated nsLTP. In the labeled-reference comparison, broad prolamin similarity strongly exceeded the best reviewed nsLTP match, and PredGPI did not support an LTPg interpretation. RsLTP20 is therefore removed from the provisional canonical set and treated as an ambiguous nsLTP/prolamin-boundary candidate. A final exclusion decision still requires stable placement in the revised positive/negative-control phylogeny and an explicit post-signal-peptide TM assessment.
"""
    (OUT / "RsLTP20_case_report_preliminary.md").write_text(report, encoding="utf-8")
    print(f"Wrote {len(output)} candidate rows to {matrix}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "build"))
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare()
    else:
        build()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
