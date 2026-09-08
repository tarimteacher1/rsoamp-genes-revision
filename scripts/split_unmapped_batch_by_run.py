#!/usr/bin/env python3
"""Split sorted category assignments and paired FASTQ by SRR run in one pass."""

from __future__ import annotations

import argparse
import subprocess
from collections import Counter
from pathlib import Path
from typing import List, Set


BUFFER_LIMIT = 4 * 1024 * 1024


class PigzWriter:
    def __init__(self, path: Path, threads: int = 1):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("wb")
        self.process = subprocess.Popen(
            ["/usr/bin/pigz", "-p", str(threads), "-c"],
            stdin=subprocess.PIPE,
            stdout=self.handle,
        )
        if self.process.stdin is None:
            raise RuntimeError(f"Could not open pigz stdin for {path}")
        self.stdin = self.process.stdin
        self.buffer = bytearray()
        self.path = path

    def write(self, data: bytes):
        self.buffer.extend(data)
        if len(self.buffer) >= BUFFER_LIMIT:
            self.flush()

    def flush(self):
        if self.buffer:
            self.stdin.write(self.buffer)
            self.buffer.clear()

    def close(self):
        self.flush()
        self.stdin.close()
        status = self.process.wait()
        self.handle.close()
        if status != 0:
            raise RuntimeError(f"pigz failed for {self.path} with status {status}")


def pigz_reader(path: Path):
    process = subprocess.Popen(["/usr/bin/pigz", "-dc", str(path)], stdout=subprocess.PIPE)
    if process.stdout is None:
        raise RuntimeError(f"Could not open pigz stdout for {path}")
    return process, process.stdout


def strip_pair_suffix(token: bytes) -> bytes:
    return token[:-2] if token.endswith((b"/1", b"/2")) else token


def run_from_id(identifier: bytes, allowed: Set[str]) -> str:
    token = identifier.split(None, 1)[0].lstrip(b"@>").decode("ascii")
    if token.endswith(("/1", "/2")):
        token = token[:-2]
    run = token.split(".", 1)[0]
    if run not in allowed:
        raise ValueError(f"Unexpected run prefix in identifier: {token}")
    return run


def split_assignments(path: Path, root: Path, runs: List[str]):
    allowed = set(runs)
    writers = {
        run: PigzWriter(root / run / "fragment_category_assignments.tsv.gz") for run in runs
    }
    counts = Counter()
    category_counts = Counter()
    previous = None
    process, stream = pigz_reader(path)
    try:
        for line in stream:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed assignment line: {line[:120]!r}")
            identifier = fields[0]
            if identifier == previous:
                raise ValueError(f"Duplicate assignment for fragment {identifier.decode('ascii')}")
            previous = identifier
            run = run_from_id(identifier, allowed)
            writers[run].write(line)
            counts[run] += 1
            category_counts[(run, fields[1].decode("ascii"))] += 1
    finally:
        stream.close()
        reader_status = process.wait()
        for writer in writers.values():
            writer.close()
    if reader_status != 0:
        raise RuntimeError(f"pigz failed while reading {path}")
    return counts, category_counts


def read_record(stream, label: str):
    header = stream.readline()
    if not header:
        return None
    sequence = stream.readline()
    plus = stream.readline()
    quality = stream.readline()
    if not sequence or not plus or not quality:
        raise ValueError(f"Truncated FASTQ record in {label}")
    return header + sequence + plus + quality, header


def split_fastq(read1: Path, read2: Path, root: Path, runs: List[str]):
    allowed = set(runs)
    writers1 = {run: PigzWriter(root / run / "residual_genome_unmapped_1.fastq.gz") for run in runs}
    writers2 = {run: PigzWriter(root / run / "residual_genome_unmapped_2.fastq.gz") for run in runs}
    counts = Counter()
    p1, s1 = pigz_reader(read1)
    p2, s2 = pigz_reader(read2)
    try:
        while True:
            rec1 = read_record(s1, str(read1))
            rec2 = read_record(s2, str(read2))
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError("R1/R2 record counts differ")
            data1, header1 = rec1
            data2, header2 = rec2
            id1 = strip_pair_suffix(header1.split(None, 1)[0].lstrip(b"@"))
            id2 = strip_pair_suffix(header2.split(None, 1)[0].lstrip(b"@"))
            if id1 != id2:
                raise ValueError(f"R1/R2 identifier mismatch: {id1!r} != {id2!r}")
            run = run_from_id(id1, allowed)
            writers1[run].write(data1)
            writers2[run].write(data2)
            counts[run] += 1
    finally:
        s1.close()
        s2.close()
        status1 = p1.wait()
        status2 = p2.wait()
        for writer in writers1.values():
            writer.close()
        for writer in writers2.values():
            writer.close()
    if status1 != 0 or status2 != 0:
        raise RuntimeError(f"pigz input failure: R1={status1}, R2={status2}")
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--counts", type=Path, required=True)
    args = parser.parse_args()

    assignment_counts, category_counts = split_assignments(
        args.assignments, args.output_root, args.runs
    )
    fastq_counts = split_fastq(args.read1, args.read2, args.output_root, args.runs)

    rows = ["run\tassignment_count\tgenome_unmapped_residual_count\tfastq_pair_count"]
    for run in args.runs:
        genome_unmapped = category_counts[(run, "genome_unmapped_residual")]
        if genome_unmapped != fastq_counts[run]:
            raise ValueError(
                f"Residual FASTQ/category mismatch for {run}: "
                f"category={genome_unmapped}, fastq={fastq_counts[run]}"
            )
        rows.append(
            f"{run}\t{assignment_counts[run]}\t{genome_unmapped}\t{fastq_counts[run]}"
        )
    args.counts.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"PASS runs={len(args.runs)} assignments={sum(assignment_counts.values())} pairs={sum(fastq_counts.values())}")


if __name__ == "__main__":
    main()
