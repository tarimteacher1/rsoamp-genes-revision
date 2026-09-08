#!/usr/bin/env python3
"""Convert TMbed three-line predictions into one auditable row per protein."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def runs(labels: str, symbols: set[str]) -> list[tuple[int, int]]:
    intervals = []
    start = None
    for index, label in enumerate(labels, start=1):
        if label in symbols and start is None:
            start = index
        elif label not in symbols and start is not None:
            intervals.append((start, index - 1))
            start = None
    if start is not None:
        intervals.append((start, len(labels)))
    return intervals


def main() -> int:
    source = Path(sys.argv[1])
    target = Path(sys.argv[2])
    lines = [line.strip() for line in source.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) % 3:
        raise ValueError(f"TMbed output does not contain complete three-line records: {source}")

    rows = []
    for offset in range(0, len(lines), 3):
        header, sequence, labels = lines[offset : offset + 3]
        if not header.startswith(">") or len(sequence) != len(labels):
            raise ValueError(f"Malformed TMbed record at line {offset + 1}")
        signal = runs(labels, {"S"})
        helices = runs(labels, {"H", "h"})
        beta = runs(labels, {"B", "b"})
        signal_end = signal[0][1] if signal and signal[0][0] <= 5 else 0
        post_signal_helices = [interval for interval in helices if interval[0] > signal_end]
        rows.append(
            {
                "protein_id": header[1:].split()[0],
                "length_aa": len(sequence),
                "signal_peptide": "YES" if signal_end else "NO",
                "signal_interval": ";".join(f"{a}-{b}" for a, b in signal),
                "tm_alpha_helix_count": len(helices),
                "tm_alpha_intervals": ";".join(f"{a}-{b}" for a, b in helices),
                "post_signal_tm_helix": "YES" if post_signal_helices else "NO",
                "post_signal_tm_intervals": ";".join(f"{a}-{b}" for a, b in post_signal_helices),
                "tm_beta_segment_count": len(beta),
                "tm_beta_intervals": ";".join(f"{a}-{b}" for a, b in beta),
                "raw_labels": labels,
            }
        )

    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} TMbed summary rows to {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
