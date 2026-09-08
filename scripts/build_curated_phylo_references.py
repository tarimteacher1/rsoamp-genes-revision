#!/usr/bin/env python3
"""Download reviewed plant family references from the official UniProt REST API."""

from __future__ import annotations

import csv
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path


REVISION = Path("/path/to/rsoamp/revision_R1_20260902")
OUT = REVISION / "02_phylogeny/reference"
OUT.mkdir(parents=True, exist_ok=True)

QUERIES = {
    "Defensin": "(reviewed:true) AND (taxonomy_id:33090) AND (xref:pfam-PF00304)",
    "Snakin_GASA": "(reviewed:true) AND (taxonomy_id:33090) AND (xref:pfam-PF02704)",
}
LENGTH_LIMITS = {"Defensin": (40, 180), "Snakin_GASA": (70, 350)}
FIELDS = "accession,id,protein_name,organism_name,organism_id,length,sequence,reviewed"


def download(query: str) -> tuple[str, str]:
    params = urllib.parse.urlencode(
        {"query": query, "format": "tsv", "fields": FIELDS, "size": 500}
    )
    url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
    request = urllib.request.Request(url, headers={"User-Agent": "Genes-RsoAMP-revision/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return url, response.read().decode("utf-8")


def main() -> int:
    manifest = []
    summary = {}
    for family, query in QUERIES.items():
        url, text = download(query)
        raw = OUT / f"{family}.uniprot_reviewed.raw.tsv"
        raw.write_text(text, encoding="utf-8")
        rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
        low, high = LENGTH_LIMITS[family]
        retained = [row for row in rows if low <= int(row["Length"]) <= high and row["Sequence"]]

        table = OUT / f"{family}.uniprot_reviewed.filtered.tsv"
        with table.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(retained)

        fasta = OUT / f"{family}.uniprot_reviewed.filtered.faa"
        with fasta.open("w", encoding="utf-8") as handle:
            for row in retained:
                organism = row["Organism (ID)"]
                handle.write(f">UNI|{row['Entry']}|taxid{organism}|{family}\n{row['Sequence']}\n")

        summary[family] = {"downloaded_reviewed": len(rows), "retained_by_length": len(retained)}
        manifest.append(
            {
                "family": family,
                "retrieval_date": date.today().isoformat(),
                "query": query,
                "url": url,
                "raw_sha256": hashlib.sha256(text.encode()).hexdigest(),
                "downloaded_reviewed": len(rows),
                "retained_by_length": len(retained),
                "length_filter_aa": f"{low}-{high}",
            }
        )

    with (OUT / "uniprot_reference_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    (OUT / "uniprot_reference_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
