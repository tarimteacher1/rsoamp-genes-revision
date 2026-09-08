#!/usr/bin/env python3
"""Assemble auditable full-precursor and cysteine-core phylogeny inputs."""

from __future__ import annotations

import csv
import hashlib
import re
from collections import defaultdict
from pathlib import Path


ORIGINAL = Path("/path/to/rsoamp")
REVISION = ORIGINAL / "revision_R1_20260902"
PHYLO = REVISION / "02_phylogeny"


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


def read_table(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def core(sequence: str) -> str:
    first = sequence.find("C")
    last = sequence.rfind("C")
    return sequence[first : last + 1] if first >= 0 and last - first + 1 >= 25 else ""


def add(rows: list[dict], family: str, taxon: str, sequence: str, source: str, source_id: str, note: str = "") -> None:
    rows.append(
        {
            "family": family,
            "taxon_id": safe(taxon),
            "source_group": source,
            "source_id": source_id,
            "length_aa": len(sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "core_length_aa": len(core(sequence)),
            "note": note,
            "sequence": sequence,
        }
    )


def load_uniprot_family(rows: list[dict], family: str) -> None:
    path = PHYLO / f"reference/{family}.uniprot_reviewed.filtered.tsv"
    seen = set()
    for row in read_table(path, "\t"):
        sequence = row["Sequence"]
        if sequence in seen:
            continue
        seen.add(sequence)
        add(
            rows,
            family,
            f"UNI_{row['Entry']}",
            sequence,
            "UniProt_reviewed_Viridiplantae",
            row["Entry"],
            f"taxid={row['Organism (ID)']}; {row['Protein names']}",
        )


def load_tamarix(rows: list[dict], family: str) -> None:
    if family == "nsLTP":
        path = REVISION / "05_comparative_genomics/predictions/tamarix_nsLTP_final_canonical.faa"
    else:
        path = REVISION / f"05_comparative_genomics/predictions/tamarix_{family}_high_confidence.faa"
    if not path.exists():
        return
    for source_id, sequence in read_fasta(path).items():
        parts = source_id.split("|")
        protein_id = parts[1] if len(parts) > 1 else source_id
        add(rows, family, f"TAU_{protein_id}", sequence, "Tamarix_austromongolica", protein_id)


def main() -> int:
    rows: list[dict] = []
    members = read_table(ORIGINAL / "tables/T1_AMP_members.csv")
    member_sequences = read_fasta(ORIGINAL / "intermediate/members_final.faa")

    # The six historical tblastn backfills are not accepted as genes before RNA/gene-model rescue.
    for member in members:
        family = member["family"]
        if family not in {"Defensin", "Snakin_GASA"} or member["source"] != "annotated":
            continue
        member_id = member["member_id"]
        add(
            rows,
            family,
            f"RSO_{member_id}",
            member_sequences[member_id],
            "Reaumuria_annotated_gene",
            member["protein_id"],
            f"member_id={member_id}",
        )

    load_uniprot_family(rows, "Defensin")
    load_uniprot_family(rows, "Snakin_GASA")
    load_tamarix(rows, "Defensin")
    load_tamarix(rows, "Snakin_GASA")

    nsltp_sequences = read_fasta(REVISION / "01_nslTP_curation/inputs/nsLTP_candidate_pool.faa")
    evidence = {
        row["protein_id"]: row
        for row in read_table(REVISION / "01_nslTP_curation/nsLTP_evidence_matrix_preliminary.tsv", "\t")
    }
    for protein_id, sequence in nsltp_sequences.items():
        row = evidence[protein_id]
        label = row["member_id"] or protein_id
        add(
            rows,
            "nsLTP",
            f"RSO_{label}",
            sequence,
            "Reaumuria_nsLTP_candidate",
            protein_id,
            row["preliminary_evidence_tier"],
        )

    reference_sets = [
        ("positive_reviewed_arabidopsis_nsLTP.faa", "NSPOS", "reviewed_nsLTP_positive"),
        ("negative_2S_seed_storage.faa", "NEG2S", "reviewed_2S_negative"),
        ("ambiguous_prolamin_umbrella.faa", "AMBPRO", "broad_prolamin_control"),
    ]
    seen_reference_sequences = set()
    for filename, prefix, group in reference_sets:
        for source_id, sequence in read_fasta(REVISION / f"01_nslTP_curation/reference/{filename}").items():
            if sequence in seen_reference_sequences:
                continue
            seen_reference_sequences.add(sequence)
            add(rows, "nsLTP", f"{prefix}_{source_id}", sequence, group, source_id)
    load_tamarix(rows, "nsLTP")

    for family in ("Defensin", "Snakin_GASA", "nsLTP"):
        family_dir = PHYLO / family
        family_dir.mkdir(parents=True, exist_ok=True)
        selected = [row for row in rows if row["family"] == family]
        taxa = [row["taxon_id"] for row in selected]
        if len(taxa) != len(set(taxa)):
            duplicates = sorted({taxon for taxon in taxa if taxa.count(taxon) > 1})
            raise ValueError(f"Duplicate taxon IDs for {family}: {duplicates}")

        for mode in ("full", "core"):
            with (family_dir / f"{family}.{mode}.faa").open("w", encoding="utf-8") as handle:
                for row in selected:
                    sequence = row["sequence"] if mode == "full" else core(row["sequence"])
                    if sequence:
                        handle.write(f">{row['taxon_id']}\n{sequence}\n")

        manifest = family_dir / f"{family}.taxon_manifest.tsv"
        columns = [key for key in selected[0] if key != "sequence"]
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows({key: row[key] for key in columns} for row in selected)
        print(f"{family}: {len(selected)} full taxa; {sum(bool(core(row['sequence'])) for row in selected)} core taxa")

    with (PHYLO / "phylogeny_inputs.sha256").open("w", encoding="utf-8") as handle:
        for path in sorted(PHYLO.glob("*/*.faa")):
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(PHYLO)}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
