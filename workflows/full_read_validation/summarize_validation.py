#!/usr/bin/env python3
"""Integrate full-read, Iso-Seq, ORF, and family evidence for R1.M6 candidates."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


ATTRIBUTE_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)"')


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"No rows and no explicit columns for {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_gtf(path: Path) -> dict[str, dict[str, object]]:
    models: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = dict(ATTRIBUTE_RE.findall(fields[8]))
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            model = models.setdefault(
                tid,
                {
                    "chrom": fields[0],
                    "strand": fields[6],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "exons": [],
                },
            )
            model["start"] = min(int(model["start"]), int(fields[3]))
            model["end"] = max(int(model["end"]), int(fields[4]))
            if fields[2] == "exon":
                model["exons"].append((int(fields[3]) - 1, int(fields[4])))
    for tid, model in models.items():
        model["exons"] = sorted(set(model["exons"]))
        if not model["exons"]:
            raise ValueError(f"No exons for transcript {tid} in {path}")
    return models


def read_fasta(path: Path) -> dict[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
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
            elif current is not None:
                records[current].append(line)
    return {key: "".join(value) for key, value in records.items()}


def transcript_from_orf(orf_id: str) -> str:
    return orf_id.rsplit("_", 1)[0]


def parse_domains(path: Path) -> dict[str, list[str]]:
    domains: dict[str, list[str]] = defaultdict(list)
    if not path.is_file():
        return domains
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.split(maxsplit=22)
            domains[fields[0]].append(fields[4])
    return domains


def parse_diamond(path: Path) -> dict[str, list[list[str]]]:
    hits: dict[str, list[list[str]]] = defaultdict(list)
    if not path.is_file():
        return hits
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) >= 14:
                hits[fields[0]].append(fields)
    for rows in hits.values():
        rows.sort(key=lambda row: -float(row[13]))
    return hits


def parse_tmbed(path: Path) -> dict[str, dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError(f"Malformed TMbed prediction file: {path}")
    results: dict[str, dict[str, object]] = {}
    for offset in range(0, len(lines), 3):
        header, sequence, labels = lines[offset : offset + 3]
        if not header.startswith(">") or len(sequence) != len(labels):
            raise ValueError(f"Malformed TMbed record near line {offset + 1}")
        identifier = header[1:].split()[0].removesuffix("|start_stop")
        signal_end = 0
        if labels.startswith("S"):
            signal_end = len(labels) - len(labels.lstrip("S"))
        post_signal_tm = any(label in {"H", "h", "B", "b"} for label in labels[signal_end:])
        results[identifier] = {
            "tmbed_signal_peptide": "YES" if signal_end else "NO",
            "tmbed_signal_end_aa": signal_end or "NA",
            "tmbed_post_signal_membrane_segment": "YES" if post_signal_tm else "NO",
            "tmbed_raw_labels": labels,
        }
    return results


def parse_predgpi(path: Path) -> dict[str, dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    records = json.loads(path.read_text(encoding="utf-8"))
    results: dict[str, dict[str, object]] = {}
    for record in records:
        identifier = str(record.get("accession", "")).split("|start_stop", 1)[0]
        if not identifier:
            continue
        features = [
            feature
            for feature in record.get("features", [])
            if feature.get("description") == "GPI-anchor"
        ]
        feature = features[0] if features else {}
        results[identifier] = {
            "predgpi_GPI_anchor": "YES" if feature else "NO",
            "predgpi_omega_site": feature.get("begin", "NA"),
            "predgpi_score": feature.get("score", "NA"),
        }
    return results


def candidate_eight_cm(sequence: str) -> tuple[str, str]:
    pattern = re.compile(r"C[^C]{1,30}C[^C]{1,35}CC[^C]{1,35}C.C[^C]{1,35}C[^C]{1,35}C")
    match = pattern.search(sequence)
    if not match:
        return "NO", ""
    positions = [index + 1 for index, aa in enumerate(match.group()) if aa == "C"]
    spacing = ",".join(str(right - left - 1) for left, right in zip(positions, positions[1:]))
    return "YES", spacing


def as_float(value: object) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def relationship_index(path: Path) -> dict[str, list[dict[str, str]]]:
    rows = read_tsv(path)
    by_ref: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("ref_id"):
            by_ref[row["ref_id"]].append(row)
    return by_ref


RELATION_RANK = {"=": 0, "c": 1, "j": 2, "o": 3, "k": 4, "m": 5, "n": 6, "x": 7, "s": 8, "i": 9, "p": 10, "u": 11}


def best_relationship(rows: list[dict[str, str]]) -> tuple[str, list[str]]:
    if not rows:
        return "NONE", []
    best_code = min((row["class_code"] for row in rows), key=lambda code: RELATION_RANK.get(code, 99))
    query_ids = sorted({row["qry_id"] for row in rows if row["class_code"] == best_code})
    return best_code, query_ids


def rna_strand(read) -> str:
    return "+" if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse) else "-"


def overlaps(blocks: list[tuple[int, int]], exons: list[tuple[int, int]]) -> bool:
    return any(left < exon_end and exon_start < right for left, right in blocks for exon_start, exon_end in exons)


def paired_support(bam, model: dict[str, object]) -> dict[str, object]:
    exons = model["exons"]
    fragments: set[str] = set()
    unique: set[str] = set()
    multiple: set[str] = set()
    covered: set[int] = set()
    introns = [(exons[index][1], exons[index + 1][0]) for index in range(len(exons) - 1)]
    junction_names: list[set[str]] = [set() for _ in introns]
    for read in bam.fetch(model["chrom"], min(a for a, _ in exons), max(b for _, b in exons)):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if rna_strand(read) != model["strand"]:
            continue
        blocks = read.get_blocks()
        if not overlaps(blocks, exons):
            continue
        fragments.add(read.query_name)
        nh = read.get_tag("NH") if read.has_tag("NH") else 0
        if nh == 1:
            unique.add(read.query_name)
            for left, right in blocks:
                for exon_start, exon_end in exons:
                    covered.update(range(max(left, exon_start), min(right, exon_end)))
            read_introns = {(blocks[index][1], blocks[index + 1][0]) for index in range(len(blocks) - 1)}
            for index, intron in enumerate(introns):
                if intron in read_introns:
                    junction_names[index].add(read.query_name)
        elif nh > 1:
            multiple.add(read.query_name)
    exon_bases = sum(end - start for start, end in exons)
    return {
        "overlapping_fragments": len(fragments),
        "unique_fragments": len(unique),
        "multimapping_fragments": len(multiple),
        "exon_bases": exon_bases,
        "unique_covered_bases": len(covered),
        "unique_exon_coverage_fraction": len(covered) / exon_bases if exon_bases else 0.0,
        "intron_count": len(introns),
        "unique_junction_fragments_by_intron": ";".join(str(len(names)) for names in junction_names) if introns else "NA",
    }


def isoseq_support(bam, model: dict[str, object]) -> dict[str, object]:
    exons = model["exons"]
    introns = [(exons[index][1], exons[index + 1][0]) for index in range(len(exons) - 1)]
    primary: set[str] = set()
    mapq20: set[str] = set()
    exact_chain: set[str] = set()
    for read in bam.fetch(model["chrom"], min(a for a, _ in exons), max(b for _, b in exons)):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if ("-" if read.is_reverse else "+") != model["strand"]:
            continue
        blocks = read.get_blocks()
        if not overlaps(blocks, exons):
            continue
        primary.add(read.query_name)
        if read.mapping_quality >= 20:
            mapq20.add(read.query_name)
            read_introns = {(blocks[index][1], blocks[index + 1][0]) for index in range(len(blocks) - 1)}
            if introns and all(intron in read_introns for intron in introns):
                exact_chain.add(read.query_name)
            elif not introns and read.reference_start <= exons[0][0] + 20 and read.reference_end >= exons[0][1] - 20:
                exact_chain.add(read.query_name)
    return {
        "isoseq_primary_overlapping_reads": len(primary),
        "isoseq_mapq20_overlapping_reads": len(mapq20),
        "isoseq_mapq20_exact_chain_or_near_full_single_exon_reads": len(exact_chain),
    }


def peptide_matches(peptide: str, related_transcripts: list[str], proteins: dict[str, str]) -> list[str]:
    related = set(related_transcripts)
    matches = []
    for oid, sequence in proteins.items():
        if transcript_from_orf(oid) in related and peptide in sequence:
            matches.append(oid)
    return sorted(matches)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-results", required=True, type=Path)
    parser.add_argument("--prep", required=True, type=Path)
    parser.add_argument("--full-results", required=True, type=Path)
    parser.add_argument("--isoseq-results", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    source = args.source_results.resolve()
    prep = args.prep.resolve()
    full = args.full_results.resolve()
    isoseq = args.isoseq_results.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    candidates = [row for row in read_tsv(source / "candidate_evidence.tsv") if row["model_origin"] == "INTERGENIC_CANDIDATE"]
    models = parse_gtf(prep / "intergenic_candidate_models.gtf")
    if {row["transcript_id"] for row in candidates} != set(models):
        raise ValueError("Candidate evidence and candidate GTF identities differ")
    subset_support = read_tsv(source / "candidate_per_sample_support.tsv")
    subset_by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in subset_support:
        subset_by_model[row["transcript_id"]].append(row)

    exon_owners: dict[tuple[str, str, int, int], set[str]] = defaultdict(set)
    for tid, model in models.items():
        for start0, end0 in model["exons"]:
            exon_owners[(str(model["chrom"]), str(model["strand"]), start0, end0)].add(tid)

    import pysam

    full_rows: list[dict[str, object]] = []
    full_by_model: dict[str, list[dict[str, object]]] = defaultdict(list)
    bam_paths = sorted((full / "align").glob("*/*.candidate_windows.bam"))
    if len(bam_paths) != 9:
        raise ValueError(f"Expected nine full-read local BAMs, observed {len(bam_paths)}")
    for bam_path in bam_paths:
        run = bam_path.parent.name
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for tid, model in models.items():
                stats = paired_support(bam, model)
                row = {"run": run, "transcript_id": tid, **stats}
                full_rows.append(row)
                full_by_model[tid].append(row)
    write_tsv(out / "candidate_full_read_per_sample_support.tsv", full_rows)

    isoseq_by_model: dict[str, dict[str, object]] = {}
    isoseq_bam = isoseq / "candidate_reads.global_remap.bam"
    if isoseq_bam.is_file() and isoseq_bam.stat().st_size:
        with pysam.AlignmentFile(str(isoseq_bam), "rb") as bam:
            for tid, model in models.items():
                isoseq_by_model[tid] = isoseq_support(bam, model)
    iso_rows = [{"transcript_id": tid, **isoseq_by_model.get(tid, {})} for tid in sorted(models)]
    write_tsv(
        out / "candidate_isoseq_read_support.tsv",
        iso_rows,
        ["transcript_id", "isoseq_primary_overlapping_reads", "isoseq_mapq20_overlapping_reads", "isoseq_mapq20_exact_chain_or_near_full_single_exon_reads"],
    )

    full_tmap = relationship_index(full / "screen/vs_subset_candidate.full_read_merged.gtf.tmap")
    iso_tmap = relationship_index(isoseq / "screen/vs_subset_candidate.isoseq_candidate_models.gtf.tmap")
    full_proteins = read_fasta(full / "screen/full_read_stop_delimited_orfs.faa")
    iso_proteins = read_fasta(isoseq / "screen/isoseq_candidate_orfs.faa")
    priority_peptides = {
        key.removesuffix("|start_stop"): value
        for key, value in read_fasta(prep / "priority_start_stop_suborfs.faa").items()
    }
    sequence_manifest = {row["orf_id"]: row for row in read_tsv(prep / "candidate_sequence_manifest.tsv")}
    tmbed = parse_tmbed(full / "screen/priority_start_stop_suborfs.tmbed.pred")
    predgpi = parse_predgpi(full / "screen/priority_start_stop_suborfs.predgpi.json")
    nsltp_similarity = {
        row["protein_id"].removesuffix("|start_stop"): row
        for row in read_tsv(full / "screen/priority_nsltp_reference_similarity.tsv")
    }
    priority_domains = parse_domains(full / "screen/priority_family_domains.domtbl")
    priority_ref_hits = parse_diamond(full / "screen/priority_vs_labeled_AMP_refs.tsv")
    proteome_hits = parse_diamond(full / "screen/priority_vs_rso_proteome.tsv")

    model_summary: list[dict[str, object]] = []
    integrated: list[dict[str, object]] = []
    for candidate in sorted(candidates, key=lambda row: (row["chrom"], int(row["start"]), row["transcript_id"])):
        tid = candidate["transcript_id"]
        oid = candidate["orf_id"]
        full_support = full_by_model[tid]
        reliable = sum(int(row["unique_fragments"]) >= 3 and float(row["unique_exon_coverage_fraction"]) >= 0.7 for row in full_support)
        subset = subset_by_model[tid]
        subset_unique = sum(int(row["unique_fragments"]) for row in subset)
        full_unique = sum(int(row["unique_fragments"]) for row in full_support)
        full_cov = [float(row["unique_exon_coverage_fraction"]) for row in full_support]
        intron_count = int(full_support[0]["intron_count"])
        exon_list = models[tid]["exons"]
        intron_lengths = [exon_list[index + 1][0] - exon_list[index][1] for index in range(len(exon_list) - 1)]
        shared_exons = sorted(
            {
                other
                for start0, end0 in exon_list
                for other in exon_owners[(str(models[tid]["chrom"]), str(models[tid]["strand"]), start0, end0)]
                if other != tid
            }
        )
        if intron_count:
            intron_sample_counts = [0] * intron_count
            for row in full_support:
                values = [int(value) for value in str(row["unique_junction_fragments_by_intron"]).split(";")]
                for index, value in enumerate(values):
                    intron_sample_counts[index] += value >= 2
            intron_gate = all(count >= 2 for count in intron_sample_counts)
        else:
            intron_sample_counts = []
            intron_gate = True
        full_class, full_query_ids = best_relationship(full_tmap.get(tid, []))
        iso_class, iso_query_ids = best_relationship(iso_tmap.get(tid, []))
        peptide = priority_peptides.get(oid, "")
        full_matches = peptide_matches(peptide, full_query_ids, full_proteins) if peptide else []
        iso_matches = peptide_matches(peptide, iso_query_ids, iso_proteins) if peptide else []
        iso_stats = isoseq_by_model.get(tid, {})

        model_row = {
            "transcript_id": tid,
            "orf_id": oid,
            "chrom": candidate["chrom"],
            "start": candidate["start"],
            "end": candidate["end"],
            "strand": candidate["strand"],
            "review_priority": candidate["review_priority"],
            "subset_total_unique_fragments": subset_unique,
            "full_read_total_unique_fragments": full_unique,
            "full_read_reliable_samples": reliable,
            "full_read_min_exon_coverage_fraction": f"{min(full_cov):.6f}",
            "full_read_mean_exon_coverage_fraction": f"{sum(full_cov) / len(full_cov):.6f}",
            "candidate_intron_count": intron_count,
            "maximum_candidate_intron_bp": max(intron_lengths) if intron_lengths else 0,
            "candidate_models_sharing_exact_exons": ";".join(shared_exons),
            "candidate_structure_review_flag": (
                "LONG_INTRON_AND_EXON_REUSE"
                if intron_lengths and max(intron_lengths) > 10000 and shared_exons
                else "EXACT_EXON_REUSE"
                if shared_exons
                else "LONG_INTRON"
                if intron_lengths and max(intron_lengths) > 10000
                else "NONE"
            ),
            "samples_with_at_least_2_unique_fragments_by_intron": ";".join(map(str, intron_sample_counts)) if intron_count else "NA",
            "full_read_intron_gate": "PASS" if intron_gate else "FAIL",
            "full_read_assembly_best_class": full_class,
            "full_read_assembly_related_transcripts": ";".join(full_query_ids),
            "exact_start_stop_peptide_in_full_read_assembly_orfs": len(full_matches),
            "isoseq_assembly_best_class": iso_class,
            "isoseq_assembly_related_transcripts": ";".join(iso_query_ids),
            "exact_start_stop_peptide_in_isoseq_assembly_orfs": len(iso_matches),
            **iso_stats,
        }
        model_summary.append(model_row)

        if candidate["review_priority"] != "FULL_READ_MODEL_VALIDATION_PRIORITY":
            continue
        features = sequence_manifest[oid]
        tm = tmbed.get(oid, {})
        gpi = predgpi.get(oid, {})
        qid = oid + "|start_stop"
        refs = priority_ref_hits.get(qid, [])
        positive = next((row for row in refs if not row[1].startswith("CONTROL_")), None)
        negative = next((row for row in refs if row[1].startswith("CONTROL_")), None)
        proteome = proteome_hits.get(qid, [])
        family_tokens = sorted(set(priority_domains.get(qid, [])))
        best_positive_score = float(positive[13]) if positive else -1.0
        best_negative_score = float(negative[13]) if negative else -1.0
        family_reference = positive[1] if positive else ""
        family_hint = (
            "Defensin" if "PF00304.27" in family_tokens or "Defensin" in family_reference or "RsDEF" in family_reference
            else "nsLTP" if "PF00234.28" in family_tokens or "nsLTP" in family_reference or "RsLTP" in family_reference
            else "UNRESOLVED"
        )
        cysteines = int(features["start_stop_cysteine_count"])
        signal_gate = tm.get("tmbed_signal_peptide") == "YES"
        model_gate = reliable >= 2 and intron_gate and full_class in {"=", "c", "j", "o"} and bool(full_matches)
        eight_cm, eight_cm_spacing = candidate_eight_cm(peptide)
        nsltp_row = nsltp_similarity.get(oid, {})
        nsltp_margin = as_float(nsltp_row.get("positive_minus_strongest_alternative_bitscore"))
        defensin_reference = bool(positive) and (
            "Defensin" in str(positive[1]) or "RsDEF" in str(positive[1])
        )
        if family_hint == "Defensin":
            family_gate = cysteines >= 6 and (
                "PF00304.27" in family_tokens
                or (defensin_reference and best_positive_score >= 30)
            )
            topology_gate = signal_gate and tm.get("tmbed_post_signal_membrane_segment") == "NO"
            family_rule_note = "Defensin-specific HMM or significant curated defensin homology plus cysteine-rich secreted precursor"
            if model_gate and topology_gate and family_gate:
                status = "SUPPORTED_POST_HOC_DEFENSIN_FAMILY_TRANSCRIPT_CANDIDATE"
            elif not model_gate:
                status = "FULL_READ_MODEL_GATE_NOT_PASSED"
            else:
                status = "DEFENSIN_FAMILY_OR_SECRETION_GATE_NOT_PASSED"
        elif family_hint == "nsLTP":
            gpi_positive = gpi.get("predgpi_GPI_anchor") == "YES"
            positive_preferred = nsltp_row.get("reference_similarity_interpretation") == "positive_reference_preferred"
            # The frozen R2 LTPg rule requires a >=50-bit positive-reference
            # margin and phylogenetic concordance. A broad PF00234 hit alone is
            # deliberately insufficient after Reviewer 1's boundary critique.
            strict_ltpg_sequence_gate = (
                signal_gate
                and eight_cm == "YES"
                and gpi_positive
                and positive_preferred
                and nsltp_margin is not None
                and nsltp_margin >= 50
            )
            family_gate = False
            topology_gate = signal_gate and (
                tm.get("tmbed_post_signal_membrane_segment") == "NO" or gpi_positive
            )
            family_rule_note = (
                "LTPg-like boundary candidate; frozen R2 sequence gate passed but candidate phylogenetic concordance is not yet established"
                if strict_ltpg_sequence_gate
                else "Broad nsLTP/prolamin candidate did not pass the frozen R1/R2 family boundary rules"
            )
            status = (
                "FULL_READ_MODEL_GATE_NOT_PASSED"
                if not model_gate
                else "NSLTP_LTPG_LIKE_BOUNDARY_CANDIDATE"
                if gpi_positive
                else "NSLTP_CORE_LIKE_BOUNDARY_CANDIDATE"
            )
        else:
            family_gate = False
            topology_gate = signal_gate and tm.get("tmbed_post_signal_membrane_segment") == "NO"
            family_rule_note = "Family identity unresolved"
            status = "FAMILY_IDENTITY_UNRESOLVED"
        integrated.append(
            {
                **model_row,
                "start_stop_length_aa": features["start_stop_length_aa"],
                "start_stop_cysteine_count": cysteines,
                "family_hint": family_hint,
                "family_HMMs_on_start_stop_peptide": ";".join(family_tokens),
                "best_labeled_AMP_reference": family_reference,
                "best_labeled_AMP_bitscore": positive[13] if positive else "NA",
                "best_boundary_control": negative[1] if negative else "",
                "best_boundary_control_bitscore": negative[13] if negative else "NA",
                "best_rso_proteome_hit": proteome[0][1] if proteome else "",
                "best_rso_proteome_pident": proteome[0][2] if proteome else "NA",
                "best_rso_proteome_query_coverage": proteome[0][10] if proteome else "NA",
                "best_rso_proteome_bitscore": proteome[0][13] if proteome else "NA",
                "candidate_8CM": eight_cm,
                "candidate_8CM_spacing": eight_cm_spacing,
                "nsltp_best_positive_reference": nsltp_row.get("best_positive_subject", ""),
                "nsltp_best_positive_bitscore": nsltp_row.get("best_positive_bitscore", "NA"),
                "nsltp_best_explicit_2S_reference": nsltp_row.get("best_explicit_2S_subject", ""),
                "nsltp_best_ambiguous_prolamin_reference": nsltp_row.get("best_ambiguous_prolamin_subject", ""),
                "nsltp_positive_minus_strongest_alternative_bitscore": nsltp_row.get("positive_minus_strongest_alternative_bitscore", "NA"),
                "nsltp_reference_similarity_interpretation": nsltp_row.get("reference_similarity_interpretation", "NOT_APPLICABLE"),
                "predgpi_GPI_anchor": gpi.get("predgpi_GPI_anchor", "NOT_RUN"),
                "predgpi_omega_site": gpi.get("predgpi_omega_site", "NA"),
                "predgpi_score": gpi.get("predgpi_score", "NA"),
                "tmbed_signal_peptide": tm.get("tmbed_signal_peptide", "NOT_RUN"),
                "tmbed_signal_end_aa": tm.get("tmbed_signal_end_aa", "NA"),
                "tmbed_post_signal_membrane_segment": tm.get("tmbed_post_signal_membrane_segment", "NOT_RUN"),
                "full_read_model_gate": "PASS" if model_gate else "FAIL",
                "family_boundary_gate": "PASS" if family_gate else "FAIL",
                "secretion_topology_gate": "PASS" if topology_gate else "FAIL",
                "family_rule_note": family_rule_note,
                "automated_evidence_status": status,
                "catalogue_action": "HOLD_FOR_MANUAL_LOCUS_AND_FAMILY_REVIEW",
            }
        )

    write_tsv(out / "candidate_full_read_model_summary.tsv", model_summary)
    write_tsv(out / "priority_candidate_integrated_evidence.tsv", integrated)
    status_counts: dict[str, int] = defaultdict(int)
    for row in integrated:
        status_counts[str(row["automated_evidence_status"])] += 1
    summary = {
        "intergenic_orf_rows_reviewed": len(candidates),
        "priority_start_stop_candidates": len(integrated),
        "full_read_local_bams": len(bam_paths),
        "automated_evidence_status_counts": dict(sorted(status_counts.items())),
        "catalogue_change": "NONE_AUTOMATIC",
        "interpretation": "Automated gates identify review priorities; they do not establish antimicrobial activity or authorize catalogue inclusion.",
    }
    (out / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
