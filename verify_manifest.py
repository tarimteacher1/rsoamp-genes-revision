"""Verify the release manifest without third-party Python dependencies."""
from pathlib import Path
import csv
import hashlib
import sys

root = Path(__file__).resolve().parent
failures = []
with (root / "MANIFEST_SHA256.tsv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
for row in rows:
    path = root / row["path"]
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
        failures.append(row["path"])
print(f"Checked {len(rows)} files; mismatches: {len(failures)}")
for path in failures:
    print(path)
sys.exit(bool(failures))
