#!/usr/bin/env python3
"""Prepare and summarize inter-species MCScanX inputs for Reaumuria and Tamarix."""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


ORIGINAL = Path("/path/to/rsoamp")
REVISION = ORIGINAL / "revision_R1_20260902"
OUT = REVISION / "05_comparative_genomics/mcscanx_rso_tamarix"
TAU = REVISION / "05_comparative_genomics/inputs/Tamarix_austromongolica"


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
    # Source proteomes encode terminal stops as either '*' or '.'.  These are
    # not residues and DIAMOND rejects '.', so normalize them without changing
    # the archived source FASTA files.
    return {
        key: "".join(value).upper().replace("*", "").replace(".", "")
        for key, value in records.items()
    }


def attrs(text: str) -> dict[str, str]:
    return dict(part.split("=", 1) for part in text.split(";") if "=" in part)


def mRNAs(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 9 and fields[2] == "mRNA":
                yield fields, attrs(fields[8])


def prepare() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    proteins = {}
    gff_rows = []
    id_rows = []

    rso_proteins = read_fasta(ORIGINAL / "data/protein.fa")
    for fields, attr in mRNAs(ORIGINAL / "data/genome.gff"):
        source_id = attr["ID"]
        if source_id not in rso_proteins:
            continue
        combined_id = f"RSO__{source_id}"
        proteins[combined_id] = rso_proteins[source_id]
        gff_rows.append((f"RSO_{fields[0]}", combined_id, int(fields[3]), int(fields[4])))
        id_rows.append((combined_id, "Reaumuria_soongarica", source_id, fields[0], fields[3], fields[4]))

    tau_proteins = read_fasta(TAU / "tau.longest_pep.fasta")
    for fields, attr in mRNAs(TAU / "tau.longest.gff3"):
        source_id = attr["ID"]
        if source_id not in tau_proteins:
            continue
        combined_id = f"TAU__{source_id}"
        proteins[combined_id] = tau_proteins[source_id]
        gff_rows.append((f"TAU_{fields[0]}", combined_id, int(fields[3]), int(fields[4])))
        id_rows.append((combined_id, "Tamarix_austromongolica", source_id, fields[0], fields[3], fields[4]))

    models = REVISION / "05_comparative_genomics/tamarix_defensin_miniprot_models.tsv"
    with models.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            combined_id = f"TAU__{row['candidate_id']}"
            proteins[combined_id] = row["sequence"]
            gff_rows.append((f"TAU_{row['chrom']}", combined_id, int(row["start"]), int(row["end"])))
            id_rows.append((combined_id, "Tamarix_austromongolica_miniprot", row["candidate_id"], row["chrom"], row["start"], row["end"]))

    with (OUT / "rso_tau.faa").open("w", encoding="utf-8") as handle:
        for protein_id, sequence in proteins.items():
            handle.write(f">{protein_id}\n{sequence}\n")
    with (OUT / "rso_tau.gff").open("w", encoding="utf-8") as handle:
        for chrom, protein_id, start, end in sorted(gff_rows):
            handle.write(f"{chrom}\t{protein_id}\t{start}\t{end}\n")
    with (OUT / "combined_id_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["combined_id", "species_source", "source_id", "chrom", "start", "end"])
        writer.writerows(id_rows)

    amp_rows = []
    final_catalogue = REVISION / "03_annotation_rescue/final_amp_primary_catalogue.tsv"
    if not final_catalogue.exists():
        raise FileNotFoundError(f"Final AMP catalogue must be frozen before MCScanX preparation: {final_catalogue}")
    with final_catalogue.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            combined_id = f"RSO__{row['protein_id']}"
            if combined_id in proteins:
                amp_rows.append((combined_id, "Reaumuria_soongarica", row["member_id"], row["family"], row["catalogue_status"]))
    with (REVISION / "05_comparative_genomics/tamarix_amp_evidence_preliminary.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["classification_status"] == "HIGH_CONFIDENCE_SECRETED_FAMILY_MEMBER" and row["final_family"] == "Snakin_GASA":
                amp_rows.append((f"TAU__{row['protein_id']}", "Tamarix_austromongolica", row["protein_id"], row["final_family"], row["classification_status"]))
    with (REVISION / "05_comparative_genomics/tamarix_nsLTP_evidence_matrix.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["final_nsLTP_classification"] == "A_CANONICAL_NSLTP":
                amp_rows.append((f"TAU__{row['protein_id']}", "Tamarix_austromongolica", row["protein_id"], "nsLTP", row["final_nsLTP_classification"]))
    with models.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            amp_rows.append((f"TAU__{row['candidate_id']}", "Tamarix_austromongolica", row["candidate_id"], "Defensin", "complete_miniprot_homology_model"))
    with (OUT / "amp_id_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["combined_id", "species", "display_id", "family", "catalogue_status"])
        writer.writerows(amp_rows)
    print(f"Prepared {len(proteins)} proteins and {len(amp_rows)} AMP entries for MCScanX")


def parse() -> None:
    amp = {}
    with (OUT / "amp_id_manifest.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            amp[row["combined_id"]] = row

    current_block = ""
    current_header = ""
    rows = []
    n_inter_blocks = 0
    with (OUT / "rso_tau.collinearity").open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("## Alignment"):
                current_header = line
                match = re.match(r"## Alignment (\d+):", line)
                current_block = match.group(1) if match else ""
                if "RSO_" in line and "TAU_" in line:
                    n_inter_blocks += 1
            elif re.match(r"^\s+\d+-\s*\d+:\s+", line) and "RSO_" in current_header and "TAU_" in current_header:
                # MCScanX separates the pair index, genes, and E value with
                # tabs.  Whitespace splitting is unstable because indices
                # below 100 are rendered as "37- 99:" but larger indices as
                # "37-100:".
                fields = [field.strip() for field in line.split("\t")]
                if len(fields) != 4:
                    raise ValueError(f"Unexpected MCScanX pair line: {line}")
                _, gene1, gene2, evalue = fields
                rso_gene = gene1 if gene1.startswith("RSO__") else gene2
                tau_gene = gene2 if gene2.startswith("TAU__") else gene1
                rso_amp = amp.get(rso_gene, {})
                tau_amp = amp.get(tau_gene, {})
                if rso_amp or tau_amp:
                    rows.append(
                        {
                            "mcscanx_block": current_block,
                            "rso_gene": rso_gene,
                            "rso_amp_id": rso_amp.get("display_id", ""),
                            "rso_amp_family": rso_amp.get("family", ""),
                            "tamarix_gene": tau_gene,
                            "tamarix_amp_id": tau_amp.get("display_id", ""),
                            "tamarix_amp_family": tau_amp.get("family", ""),
                            "pair_evalue": evalue,
                            "direct_amp_to_amp_anchor": "YES" if rso_amp and tau_amp else "NO",
                        }
                    )

    target = REVISION / "05_comparative_genomics/rso_tamarix_amp_synteny.tsv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        columns = [
            "mcscanx_block", "rso_gene", "rso_amp_id", "rso_amp_family", "tamarix_gene",
            "tamarix_amp_id", "tamarix_amp_family", "pair_evalue", "direct_amp_to_amp_anchor",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    direct = sum(row["direct_amp_to_amp_anchor"] == "YES" for row in rows)
    print(f"Parsed {n_inter_blocks} inter-species collinear blocks; {len(rows)} AMP-involving anchors; {direct} direct AMP-to-AMP anchors")


def main() -> int:
    if sys.argv[1] == "prepare":
        prepare()
    elif sys.argv[1] == "parse":
        parse()
    else:
        raise ValueError(sys.argv[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
