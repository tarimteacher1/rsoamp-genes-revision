#!/usr/bin/env python3
"""Summarize family evidence for M6 models overlapping existing annotation."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current: str | None = None
    with path.open(encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                records[current].append(line)
    return {identifier: "".join(parts) for identifier, parts in records.items()}


def parse_domains(path: Path) -> dict[str, list[tuple[str, str]]]:
    domains: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.split(maxsplit=22)
            domains[fields[0]].append((fields[3], fields[4]))
    return domains


def parse_hits(path: Path) -> dict[str, list[list[str]]]:
    hits: dict[str, list[list[str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 14:
                hits[fields[0]].append(fields)
    for rows in hits.values():
        rows.sort(key=lambda row: -float(row[13]))
    return hits


def parse_tmbed(path: Path) -> dict[str, dict[str, object]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError(f"Malformed TMbed file: {path}")
    result: dict[str, dict[str, object]] = {}
    for offset in range(0, len(lines), 3):
        header, sequence, labels = lines[offset : offset + 3]
        if not header.startswith(">") or len(sequence) != len(labels):
            raise ValueError(f"Malformed TMbed record near line {offset + 1}")
        identifier = header[1:].split()[0].removesuffix("|start_stop")
        signal_end = len(labels) - len(labels.lstrip("S")) if labels.startswith("S") else 0
        post_signal_tm = any(label in {"H", "h", "B", "b"} for label in labels[signal_end:])
        result[identifier] = {
            "tmbed_signal_peptide": "YES" if signal_end else "NO",
            "tmbed_signal_end_aa": signal_end or "NA",
            "tmbed_post_signal_membrane_segment": "YES" if post_signal_tm else "NO",
        }
    return result


def parse_predgpi(path: Path) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for record in json.loads(path.read_text(encoding="utf-8")):
        identifier = str(record.get("accession", "")).split("|start_stop", 1)[0]
        anchors = [feature for feature in record.get("features", []) if feature.get("description") == "GPI-anchor"]
        feature = anchors[0] if anchors else {}
        result[identifier] = {
            "predgpi_GPI_anchor": "YES" if feature else "NO",
            "predgpi_omega_site": feature.get("begin", "NA"),
            "predgpi_score": feature.get("score", "NA"),
        }
    return result


def eight_cm(sequence: str) -> tuple[str, str]:
    pattern = re.compile(r"C[^C]{1,30}C[^C]{1,35}CC[^C]{1,35}C.C[^C]{1,35}C[^C]{1,35}C")
    match = pattern.search(sequence)
    if not match:
        return "NO", ""
    positions = [index + 1 for index, aa in enumerate(match.group()) if aa == "C"]
    return "YES", ",".join(str(right - left - 1) for left, right in zip(positions, positions[1:]))


def as_float(value: object) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def family_hint(domain_names: set[str], best_reference: str) -> str:
    if "Gamma-thionin" in domain_names or "Defensin" in best_reference or "RsDEF" in best_reference:
        return "Defensin"
    if "Snakin" in domain_names or "GASA" in domain_names or "Snakin" in best_reference or "RsSNA" in best_reference:
        return "Snakin_GASA"
    if "Tryp_alpha_amyl" in domain_names or "nsLTP" in best_reference or "RsLTP" in best_reference:
        return "nsLTP_prolamin_boundary"
    if "Chitin_bind_1" in domain_names:
        return "chitin_binding_long_protein"
    return "UNRESOLVED"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prep", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    prep = args.prep.resolve()
    results = args.results.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(prep / "overlap_model_manifest.tsv")
    sequences = {
        identifier.removesuffix("|start_stop"): sequence
        for identifier, sequence in read_fasta(prep / "overlap_complete_start_stop_suborfs.faa").items()
    }
    domains = parse_domains(results / "overlap_family_domains.domtbl")
    references = parse_hits(results / "overlap_vs_labeled_AMP_refs.tsv")
    proteome = parse_hits(results / "overlap_vs_rso_proteome.tsv")
    nsltp = {
        row["protein_id"].removesuffix("|start_stop"): row
        for row in read_tsv(results / "overlap_nsltp_reference_similarity.tsv")
    }
    tmbed = parse_tmbed(results / "overlap_start_stop_suborfs.tmbed.pred")
    predgpi = parse_predgpi(results / "overlap_start_stop_suborfs.predgpi.json")

    reviewed: list[dict[str, object]] = []
    for row in manifest:
        oid = row["orf_id"]
        qid = oid + "|start_stop"
        sequence = sequences.get(oid, "")
        domain_pairs = domains.get(qid, [])
        domain_names = {name for name, _ in domain_pairs}
        ref_hits = references.get(qid, [])
        best_ref = ref_hits[0] if ref_hits else None
        best_ref_id = best_ref[1] if best_ref else ""
        family = family_hint(domain_names, best_ref_id)
        rso_hits = proteome.get(qid, [])
        best_rso = rso_hits[0] if rso_hits else None
        tm = tmbed.get(oid, {})
        gpi = predgpi.get(oid, {})
        nsl = nsltp.get(oid, {})
        margin = as_float(nsl.get("positive_minus_strongest_alternative_bitscore"))
        cm, spacing = eight_cm(sequence) if sequence else ("NOT_ASSESSED", "")
        complete = row["complete_start_stop_suborf"] == "YES"
        exact_reference_peptide = bool(
            best_rso
            and best_rso[1] == row["reference_id"]
            and float(best_rso[2]) == 100.0
            and float(best_rso[10]) == 100.0
            and int(best_rso[4]) == int(best_rso[5])
        )
        r2_sequence_gate = bool(
            complete
            and family == "nsLTP_prolamin_boundary"
            and tm.get("tmbed_signal_peptide") == "YES"
            and cm == "YES"
            and gpi.get("predgpi_GPI_anchor") == "YES"
            and nsl.get("reference_similarity_interpretation") == "positive_reference_preferred"
            and margin is not None
            and margin >= 50
        )

        frozen_class = row["frozen_nsltp_final_classification"]
        if not complete:
            status = "INCOMPLETE_OR_TOO_SHORT_ORF"
            action = "DO_NOT_ADD"
            reason = "No complete start-to-stop ORF meeting the discovery length gate"
        elif frozen_class:
            status = "REDISCOVERED_FROZEN_NSLTP_BOUNDARY_OR_EXCLUSION"
            action = "RETAIN_FROZEN_CLASSIFICATION"
            reason = "The model maps to a sequence already audited under the frozen nsLTP rules; read support does not override family-boundary evidence"
        elif family == "chitin_binding_long_protein":
            status = "NON_TARGET_LONG_CHITIN_BINDING_PROTEIN"
            action = "DO_NOT_ADD"
            reason = "A long chitin-binding-domain protein is not a small hevein-like precursor"
        elif family == "Defensin":
            if sequence.count("C") < 6 or "Gamma-thionin" not in domain_names or tm.get("tmbed_signal_peptide") != "YES":
                status = "WEAK_DEFENSIN_SIMILARITY_WITHOUT_FAMILY_GATE"
                action = "DO_NOT_ADD"
                reason = "Weak local homology lacks the required cysteine/domain/secretion convergence"
            else:
                status = "DEFENSIN_FAMILY_SEQUENCE_CANDIDATE_REQUIRES_LOCUS_REVIEW"
                action = "HOLD_FOR_MANUAL_REVIEW"
                reason = "Sequence-level defensin gates pass, but the model overlaps existing annotation and is not an independent locus call"
        elif family == "Snakin_GASA":
            if sequence.count("C") < 6 or tm.get("tmbed_signal_peptide") != "YES":
                status = "WEAK_SNAKIN_SIMILARITY_WITHOUT_FAMILY_GATE"
                action = "DO_NOT_ADD"
                reason = "Weak local homology lacks a cysteine-rich secreted Snakin/GASA precursor"
            else:
                status = "SNAKIN_FAMILY_SEQUENCE_CANDIDATE_REQUIRES_LOCUS_REVIEW"
                action = "HOLD_FOR_MANUAL_REVIEW"
                reason = "Sequence-level Snakin/GASA evidence requires independent locus review"
        elif family == "nsLTP_prolamin_boundary":
            if r2_sequence_gate:
                status = "NSLTP_R2_SEQUENCE_GATE_ONLY_WITHOUT_PHYLOGENETIC_OVERRIDE"
                action = "HOLD_OUTSIDE_PRIMARY_CATALOGUE"
                reason = "The R2 sequence gate passes, but no new-locus or phylogenetic evidence authorizes overriding the frozen boundary decision"
            else:
                status = "BROAD_NSLTP_PROLAMIN_EVIDENCE_ONLY"
                action = "DO_NOT_ADD"
                reason = "Broad PF00234/prolamin compatibility does not meet the frozen canonical nsLTP rules"
        else:
            status = "NO_TARGET_FAMILY_GATE"
            action = "DO_NOT_ADD"
            reason = "No convergent target-family evidence"

        reviewed.append(
            {
                **row,
                "best_rso_proteome_hit": best_rso[1] if best_rso else "",
                "best_rso_proteome_pident": best_rso[2] if best_rso else "NA",
                "best_rso_proteome_query_coverage": best_rso[10] if best_rso else "NA",
                "exact_reference_peptide": "YES" if exact_reference_peptide else "NO",
                "family_hint": family,
                "family_HMMs": ";".join(f"{name}|{accession}" for name, accession in domain_pairs),
                "best_labeled_AMP_reference": best_ref_id,
                "best_labeled_AMP_bitscore": best_ref[13] if best_ref else "NA",
                "candidate_8CM": cm,
                "candidate_8CM_spacing": spacing,
                "tmbed_signal_peptide": tm.get("tmbed_signal_peptide", "NOT_ASSESSED"),
                "tmbed_signal_end_aa": tm.get("tmbed_signal_end_aa", "NA"),
                "tmbed_post_signal_membrane_segment": tm.get("tmbed_post_signal_membrane_segment", "NOT_ASSESSED"),
                "predgpi_GPI_anchor": gpi.get("predgpi_GPI_anchor", "NOT_ASSESSED"),
                "predgpi_omega_site": gpi.get("predgpi_omega_site", "NA"),
                "nsltp_positive_minus_alternative_bitscore": nsl.get("positive_minus_strongest_alternative_bitscore", "NA"),
                "nsltp_reference_similarity_interpretation": nsl.get("reference_similarity_interpretation", "NOT_ASSESSED"),
                "frozen_R2_sequence_gate": "PASS" if r2_sequence_gate else "FAIL" if complete and family == "nsLTP_prolamin_boundary" else "NOT_APPLICABLE",
                "family_review_status": status,
                "final_catalogue_action": action,
                "decision_reason": reason,
            }
        )

    output = out / "overlap_model_family_review.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(reviewed[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(reviewed)

    status_counts = Counter(str(row["family_review_status"]) for row in reviewed)
    summary = {
        "models_reviewed": len(reviewed),
        "complete_start_stop_ORFs": sum(row["complete_start_stop_suborf"] == "YES" for row in reviewed),
        "family_review_status_counts": dict(sorted(status_counts.items())),
        "automatic_primary_catalogue_additions": 0,
        "interpretation": "Overlapping/altered models were audited against frozen family rules; transcript support alone did not create a new locus or establish antimicrobial activity.",
    }
    (out / "overlap_model_family_review_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
