#!/usr/bin/env python3
"""Convert paired gzipped FASTQ files to one mate-labeled gzipped FASTA."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--read1", required=True, type=Path)
    parser.add_argument("--read2", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def convert(path: Path, mate: int, output) -> int:
    count = 0
    with gzip.open(path, "rt", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline()
            if not plus or not quality or not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ: {path}")
            read_id = header[1:].split()[0]
            if read_id.endswith(("/1", "/2")):
                read_id = read_id[:-2]
            output.write(f">{read_id}/{mate}\n{sequence}\n")
            count += 1
    return count


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt", encoding="ascii", compresslevel=4) as output:
        n1 = convert(args.read1, 1, output)
        n2 = convert(args.read2, 2, output)
    if n1 != n2:
        raise ValueError(f"Paired FASTQ count mismatch: read1={n1} read2={n2}")
    print(f"PASS FASTA fragments={n1} reads={n1 + n2}")


if __name__ == "__main__":
    main()
