#!/usr/bin/env python3
"""Summarize Salmon input processing and mapping gates from meta_info.json."""

from __future__ import annotations

import csv
import json
from pathlib import Path


REVISION = Path("/path/to/rsoamp/revision_R1_20260902")
OUT = REVISION / "06_expression_reanalysis"


def main() -> int:
    rows = []
    with (OUT / "rnaseq_sample_sheet.tsv").open(encoding="utf-8") as handle:
        for sample in csv.DictReader(handle, delimiter="\t"):
            run = sample["run"]
            meta = json.loads((OUT / f"salmon_quant/{run}/aux_info/meta_info.json").read_text())
            library_types = meta.get("library_types", [])
            if isinstance(library_types, str):
                library_types = [library_types]
            rows.append(
                {
                    **sample,
                    "num_processed": meta.get("num_processed"),
                    "num_mapped": meta.get("num_mapped"),
                    "percent_mapped": meta.get("percent_mapped"),
                    "library_types": ";".join(library_types),
                    "mapping_gate": "PASS" if meta.get("num_processed", 0) >= 27_000_000 else "FAIL",
                }
            )
    with (OUT / "salmon_quant_integrity.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if any(row["mapping_gate"] != "PASS" for row in rows):
        raise SystemExit("At least one Salmon run processed fewer than 27 million fragments")
    print(f"PASS Salmon integrity for {len(rows)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
