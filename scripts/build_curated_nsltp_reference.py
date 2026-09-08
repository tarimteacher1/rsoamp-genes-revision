#!/usr/bin/env python3
"""Build explicit positive and negative nsLTP reference sets for revision R1."""

from __future__ import annotations

import csv
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path


ORIGINAL = Path("/path/to/rsoamp")
REVISION = ORIGINAL / "revision_R1_20260902"
OUT = REVISION / "01_nslTP_curation/reference"
UNIPROT_QUERY = '(organism_id:3702) AND (reviewed:true) AND (protein_name:"Non-specific lipid-transfer protein")'
UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/stream"
EXPLICIT_NEGATIVE = ("2s seed storage protein", "2s albumin seed storage protein")
BROAD_AMBIGUOUS = (
    "bifunctional inhibitor/plant lipid transfer protein/seed storage",
    "bifunctional inhibitor/lipid-transfer protein/seed storage",
)


def get(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Genes-R1-reproducibility-audit/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode("utf-8")


def read_fasta_text(text: str) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, list[str]]] = {}
    current = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:]
            parts = header.split("|", 2)
            accession = parts[1] if len(parts) >= 3 else header.split()[0]
            current = accession
            records[current] = (header, [])
        elif current:
            records[current][1].append(line.replace("*", ""))
    return {accession: (header, "".join(seq)) for accession, (header, seq) in records.items()}


def write_fasta(path: Path, records: list[dict], prefix: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in sorted(records, key=lambda row: row["accession"]):
            handle.write(f">{prefix}|{record['accession']} {record['description']}\n")
            sequence = record["sequence"]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    encoded = urllib.parse.urlencode({"query": UNIPROT_QUERY, "format": "fasta"})
    positive_url = f"{UNIPROT_BASE}?{encoded}"
    positive_text = get(positive_url)
    positives_raw = read_fasta_text(positive_text)

    positive_records = []
    for accession, (header, sequence) in positives_raw.items():
        description = header.split(" ", 1)[1] if " " in header else header
        category = "positive_reviewed_LTPg" if "GPI-anchored" in description else "positive_reviewed_nsLTP"
        positive_records.append(
            {
                "accession": accession,
                "description": description,
                "sequence": sequence,
                "evidence_category": category,
                "source": "UniProtKB reviewed Arabidopsis query",
                "source_url": positive_url,
            }
        )

    old_records = read_fasta_text((ORIGINAL / "data/ref/nsLTP_PF00234.fa").read_text(encoding="utf-8"))
    negative_records = []
    ambiguous_records = []
    old_audit = []
    for accession, (header, sequence) in old_records.items():
        description = header.split(" ", 1)[1] if " " in header else header
        lowered = description.lower()
        if any(term in lowered for term in EXPLICIT_NEGATIVE):
            category = "negative_2S_seed_storage"
            bucket = negative_records
        elif any(term in lowered for term in BROAD_AMBIGUOUS):
            category = "ambiguous_prolamin_umbrella"
            bucket = ambiguous_records
        elif "non-specific lipid-transfer protein" in lowered:
            category = "historical_named_nsLTP_not_used_as_curated_positive"
            bucket = None
        else:
            category = "historical_unresolved_reference"
            bucket = None
        record = {
            "accession": accession,
            "description": description,
            "sequence": sequence,
            "evidence_category": category,
            "source": "historical PF00234 reference FASTA",
            "source_url": str(ORIGINAL / "data/ref/nsLTP_PF00234.fa"),
        }
        old_audit.append(record)
        if bucket is not None:
            bucket.append(record)

    q29 = json.loads(get("https://rest.uniprot.org/uniprotkb/Q29QA0.json"))
    q29_record = {
        "accession": q29["primaryAccession"],
        "description": q29["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
        + " | " + q29["proteinDescription"]["submissionNames"][1]["fullName"]["value"],
        "sequence": q29["sequence"]["value"],
        "evidence_category": "ambiguous_prolamin_umbrella_old_RsLTP20_best_hit",
        "source": "UniProtKB Q29QA0 live record",
        "source_url": "https://rest.uniprot.org/uniprotkb/Q29QA0.json",
    }
    if all(record["sequence"] != q29_record["sequence"] for record in ambiguous_records):
        ambiguous_records.append(q29_record)

    write_fasta(OUT / "positive_reviewed_arabidopsis_nsLTP.faa", positive_records, "POS")
    write_fasta(OUT / "negative_2S_seed_storage.faa", negative_records, "NEG2S")
    write_fasta(OUT / "ambiguous_prolamin_umbrella.faa", ambiguous_records, "AMBPRO")
    combined = positive_records + negative_records + ambiguous_records
    prefixes = {
        **{record["accession"]: "POS" for record in positive_records},
        **{record["accession"]: "NEG2S" for record in negative_records},
        **{record["accession"]: "AMBPRO" for record in ambiguous_records},
    }
    with (OUT / "combined_labeled_reference.faa").open("w", encoding="utf-8") as handle:
        for record in sorted(combined, key=lambda row: (prefixes[row["accession"]], row["accession"])):
            handle.write(f">{prefixes[record['accession']]}|{record['accession']} {record['description']}\n{record['sequence']}\n")

    manifest_fields = ["accession", "description", "length_aa", "evidence_category", "source", "source_url", "sequence_sha256"]
    with (OUT / "reference_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
        writer.writeheader()
        for record in sorted(combined, key=lambda row: (row["evidence_category"], row["accession"])):
            writer.writerow(
                {
                    **{key: record[key] for key in ("accession", "description", "evidence_category", "source", "source_url")},
                    "length_aa": len(record["sequence"]),
                    "sequence_sha256": hashlib.sha256(record["sequence"].encode()).hexdigest(),
                }
            )

    with (OUT / "historical_reference_contamination_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "description", "evidence_category"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(sorted(old_audit, key=lambda row: row["accession"]))

    summary = {
        "positive_reviewed": len(positive_records),
        "negative_2S": len(negative_records),
        "ambiguous_prolamin": len(ambiguous_records),
        "historical_total": len(old_audit),
        "historical_explicit_negative_or_ambiguous": sum(
            record["evidence_category"] in {"negative_2S_seed_storage", "ambiguous_prolamin_umbrella"}
            for record in old_audit
        ),
        "old_RsLTP20_best_hit_Q29QA0_category": q29_record["evidence_category"],
    }
    (OUT / "reference_build_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
