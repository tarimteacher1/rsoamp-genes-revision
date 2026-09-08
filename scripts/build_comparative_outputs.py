#!/usr/bin/env python3
"""Summarize same-rule AMP counts, reciprocal hits and RSO-Tamarix synteny."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


BASE = Path("/path/to/rsoamp/revision_R1_20260902")
DIR = BASE / "05_comparative_genomics"
MCSCAN = DIR / "mcscanx_rso_tamarix"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rso = read_tsv(BASE / "03_annotation_rescue/final_amp_primary_catalogue.tsv")
    tau_ns = [
        row for row in read_tsv(DIR / "tamarix_nsLTP_evidence_matrix.tsv")
        if row["final_nsLTP_classification"] == "A_CANONICAL_NSLTP"
    ]
    tau_prelim = read_tsv(DIR / "tamarix_amp_evidence_preliminary.tsv")
    tau_sna = [
        row for row in tau_prelim
        if row["final_family"] == "Snakin_GASA" and row["classification_status"] == "HIGH_CONFIDENCE_SECRETED_FAMILY_MEMBER"
    ]
    tau_def = read_tsv(DIR / "tamarix_defensin_miniprot_models.tsv")
    if len(tau_def) != 5 or any(float(row["query_coverage"]) < 0.95 for row in tau_def):
        raise ValueError("Tamarix defensin rescue table does not contain five near-full-length homology models")

    rso_amp = {
        f"RSO__{row['protein_id']}": {"display_id": row["member_id"], "family": row["family"]}
        for row in rso if row["historical_source"] == "annotated"
    }
    tau_amp = {
        **{f"TAU__{row['protein_id']}": {"display_id": row["protein_id"], "family": "nsLTP"} for row in tau_ns},
        **{f"TAU__{row['protein_id']}": {"display_id": row["protein_id"], "family": "Snakin_GASA"} for row in tau_sna},
        **{f"TAU__{row['candidate_id']}": {"display_id": row["candidate_id"], "family": "Defensin"} for row in tau_def},
    }

    rso_counts = Counter(row["family"] for row in rso)
    tau_counts = Counter(item["family"] for item in tau_amp.values())
    family_rows = []
    normalized_rows = []
    for family in ("Defensin", "Snakin_GASA", "nsLTP"):
        family_rows.extend(
            [
                {"species": "Reaumuria soongarica", "family": family, "high_confidence_count": rso_counts[family], "ambiguous_not_counted": 0, "protein_denominator": 21791},
                {"species": "Tamarix austromongolica", "family": family, "high_confidence_count": tau_counts[family], "ambiguous_not_counted": 43 - len(tau_ns) if family == "nsLTP" else 0, "protein_denominator": 22374},
            ]
        )
        normalized_rows.extend(
            [
                {"species": "Reaumuria soongarica", "family": family, "high_confidence_count": rso_counts[family], "per_10000_predicted_proteins": f"{rso_counts[family] / 21791 * 10000:.4f}"},
                {"species": "Tamarix austromongolica", "family": family, "high_confidence_count": tau_counts[family], "per_10000_predicted_proteins": f"{tau_counts[family] / 22374 * 10000:.4f}"},
            ]
        )
    write_tsv(DIR / "family_count_comparison.tsv", family_rows)
    write_tsv(DIR / "normalized_family_counts.tsv", normalized_rows)

    best_cross = {}
    with (MCSCAN / "rso_tau.blast").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query, subject = fields[0], fields[1]
            if query.split("__", 1)[0] == subject.split("__", 1)[0]:
                continue
            bitscore = float(fields[11])
            evalue = float(fields[10])
            identity = float(fields[2])
            aligned_length = int(fields[3])
            current = best_cross.get(query)
            if current is None or (bitscore, -evalue) > (current["bitscore"], -current["evalue"]):
                best_cross[query] = {
                    "subject": subject, "bitscore": bitscore, "evalue": evalue,
                    "pident": identity, "alignment_length": aligned_length,
                }

    all_amp = {**rso_amp, **tau_amp}
    orthogroups = []
    seen_pairs = set()
    for query, query_info in sorted(all_amp.items()):
        hit = best_cross.get(query)
        if not hit:
            continue
        subject = hit["subject"]
        reverse = best_cross.get(subject)
        if not reverse or reverse["subject"] != query:
            continue
        pair = tuple(sorted((query, subject)))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        rso_id = query if query.startswith("RSO__") else subject
        tau_id = subject if query.startswith("RSO__") else query
        tau_info = tau_amp.get(tau_id, {})
        rso_info = rso_amp.get(rso_id, {})
        orthogroups.append(
            {
                "orthogroup_id": f"RBH_{len(orthogroups) + 1:04d}",
                "rso_combined_id": rso_id,
                "rso_display_id": rso_info.get("display_id", rso_id.replace("RSO__", "")),
                "rso_is_final_AMP": "YES" if rso_id in rso_amp else "NO",
                "rso_family": rso_info.get("family", ""),
                "tamarix_combined_id": tau_id,
                "tamarix_display_id": tau_info.get("display_id", tau_id.replace("TAU__", "")),
                "tamarix_is_final_AMP": "YES" if tau_id in tau_amp else "NO",
                "tamarix_family": tau_info.get("family", ""),
                "amp_family_concordant": "YES" if rso_info and tau_info and rso_info["family"] == tau_info["family"] else "NO" if rso_info and tau_info else "NOT_BOTH_AMP",
                "pident": hit["pident"],
                "alignment_length": hit["alignment_length"],
                "evalue": hit["evalue"],
                "bitscore": hit["bitscore"],
                "orthology_definition": "reciprocal best cross-species DIAMOND hit; pairwise screening orthogroup",
            }
        )
    write_tsv(DIR / "orthogroups.tsv", orthogroups)

    synteny = read_tsv(DIR / "rso_tamarix_amp_synteny.tsv")
    write_tsv(DIR / "synteny_pairs.tsv", synteny)
    syntenic_amp_ids = set()
    for row in synteny:
        if row["rso_amp_id"]:
            syntenic_amp_ids.add(f"RSO_DISPLAY::{row['rso_amp_id']}")
        if row["tamarix_amp_id"]:
            syntenic_amp_ids.add(f"TAU_DISPLAY::{row['tamarix_amp_id']}")
    rbh_amp_ids = set()
    for row in orthogroups:
        if row["rso_is_final_AMP"] == "YES":
            rbh_amp_ids.add(f"RSO_DISPLAY::{row['rso_display_id']}")
        if row["tamarix_is_final_AMP"] == "YES":
            rbh_amp_ids.add(f"TAU_DISPLAY::{row['tamarix_display_id']}")

    candidate_rows = []
    for species_prefix, amp_map in (("RSO", rso_amp), ("TAU", tau_amp)):
        for combined_id, info in sorted(amp_map.items()):
            key = f"{species_prefix}_DISPLAY::{info['display_id']}"
            no_rbh = key not in rbh_amp_ids
            no_synteny = key not in syntenic_amp_ids
            candidate_rows.append(
                {
                    "species": "Reaumuria soongarica" if species_prefix == "RSO" else "Tamarix austromongolica",
                    "combined_id": combined_id,
                    "display_id": info["display_id"],
                    "family": info["family"],
                    "reciprocal_best_hit_detected": "NO" if no_rbh else "YES",
                    "amp_in_mcscanx_anchor_detected": "NO" if no_synteny else "YES",
                    "screening_status": "NO_RBH_OR_SYNTENY_CANDIDATE" if no_rbh and no_synteny else "COMPARATIVE_SUPPORT_PRESENT",
                    "interpretation_limit": "Absence of pairwise RBH/synteny does not establish lineage specificity or gene-family loss",
                }
            )
    write_tsv(DIR / "lineage_specific_candidates.tsv", candidate_rows)

    metadata_rows = [
        {
            "species": "Reaumuria soongarica", "family": "Tamaricaceae", "order": "Caryophyllales",
            "assembly_accession": "JBEBFM000000000", "assembly_level": "Chromosome-level",
            "assembly_size_bp": 1281161807, "predicted_proteins": 21791, "busco_complete_percent": 97.5,
            "data_source": "10.6084/m9.figshare.25533064.v2", "publication_doi": "10.1038/s41597-024-03644-y",
            "comparison_role": "focal species",
        },
        {
            "species": "Tamarix austromongolica", "family": "Tamaricaceae", "order": "Caryophyllales",
            "assembly_accession": "GCA_039764185.1", "assembly_level": "Chromosome-level",
            "assembly_size_bp": 1327523198, "predicted_proteins": 22374, "busco_complete_percent": 98.2,
            "data_source": "10.6084/m9.figshare.25106726.v1", "publication_doi": "10.1093/dnares/dsae021",
            "comparison_role": "same-family close comparator",
        },
    ]
    write_tsv(DIR / "species_metadata.tsv", metadata_rows)

    direct_amp = sum(row["direct_amp_to_amp_anchor"] == "YES" for row in synteny)
    concordant_rbh = sum(row["amp_family_concordant"] == "YES" for row in orthogroups)
    report = f"""# Close-species comparative genomics summary

Status: PASS

Tamarix austromongolica was selected as the close comparator because it is a chromosome-level genome from the same family (Tamaricaceae), with a similarly sized predicted proteome and public sequence/GFF provenance. The same HMM, InterPro, secretion, cysteine-framework, PredGPI, TMbed, and labeled-reference gates were applied to nsLTP candidates before counts were compared.

Observed high-confidence counts were R. soongarica: Defensin={rso_counts['Defensin']}, Snakin/GASA={rso_counts['Snakin_GASA']}, nsLTP={rso_counts['nsLTP']}; and T. austromongolica: Defensin={tau_counts['Defensin']}, Snakin/GASA={tau_counts['Snakin_GASA']}, nsLTP={tau_counts['nsLTP']}. Counts are also reported per 10,000 predicted proteins.

Pairwise DIAMOND screening identified {len(orthogroups)} reciprocal-best-hit groups involving at least one final AMP, including {concordant_rbh} groups in which both partners were final AMPs of the same family. MCScanX yielded {len(synteny)} AMP-involving interspecies anchors, of which {direct_amp} directly paired final AMP candidates in both species.

These results support discussion of higher or lower *observed* family representation and conserved candidates. They do not by themselves prove lineage-specific expansion or loss; such claims require concordant phylogeny, duplication and synteny evidence and remain sensitive to the five Tamarix defensin homology models that were absent from the official annotation.
"""
    (DIR / "comparative_summary.md").write_text(report, encoding="utf-8")
    selection = """# Species selection report

Tamarix austromongolica was chosen before examining AMP family counts. It is a same-family Tamaricaceae comparator with a chromosome-level assembly (GCA_039764185.1), 22,374 longest-protein models, public GFF/proteome/genome files, checksums, and a primary genome report (doi:10.1093/dnares/dsae021). Arabidopsis thaliana is retained only as a distant reference for continuity with the submitted analysis.

The comparator is suitable for family composition, reciprocal homology and chromosome-scale synteny. Annotation sensitivity is reported explicitly: five complete defensin homology models were recovered from intergenic genome sequence and are not represented as official annotated genes.
"""
    (DIR / "species_selection_report.md").write_text(selection, encoding="utf-8")
    (DIR / "COMPARATIVE_GENOMICS_COMPLETE.PASS").write_text(
        f"PASS\trbh_groups={len(orthogroups)}\tamp_synteny_rows={len(synteny)}\tdirect_amp_anchors={direct_amp}\n",
        encoding="utf-8",
    )
    print(report)


if __name__ == "__main__":
    main()
