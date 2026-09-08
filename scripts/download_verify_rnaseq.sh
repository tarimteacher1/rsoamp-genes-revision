#!/usr/bin/env bash
# Download the nine public salt-treatment RNA-seq runs and enforce input integrity.
set -euo pipefail

REV="${REV:-/path/to/rsoamp/revision_R1_20260902}"
AUDIT="$REV/00_audit"
OUT="$REV/04_unmapped_reads/raw_fastq"
LOGDIR="$REV/04_unmapped_reads/logs"
ENA="$AUDIT/ena_run_metadata_20260902.tsv"
SEQKIT="${SEQKIT:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/seqkit}"
JOBS="${JOBS:-3}"

mkdir -p "$OUT" "$LOGDIR"
test -s "$ENA"
test -x "$SEQKIT"

MANIFEST="$AUDIT/rna_fastq_download_manifest.tsv"
STATUS="$AUDIT/rna_fastq_download_status.tsv"
DOWNLOAD_LOG="$LOGDIR/rna_fastq_download.log"
: > "$DOWNLOAD_LOG"

python3 - "$ENA" "$OUT" > "$MANIFEST" <<'PY'
import csv
import sys
from pathlib import Path

ena, out = Path(sys.argv[1]), Path(sys.argv[2])
wanted = {
    "SRR27540875", "SRR27540876", "SRR27540877",
    "SRR27540878", "SRR27540879", "SRR27540880",
    "SRR27540881", "SRR27540883", "SRR27540884",
}
print("run\tmate\turl\texpected_md5\texpected_bytes\toutput")
with ena.open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        run = row["run_accession"]
        if run not in wanted:
            continue
        urls = row["fastq_ftp"].split(";")
        md5s = row["fastq_md5"].split(";")
        sizes = row["fastq_bytes"].split(";")
        if not (len(urls) == len(md5s) == len(sizes) == 2):
            raise SystemExit(f"Expected two mates for {run}")
        for mate, (url, md5, size) in enumerate(zip(urls, md5s, sizes), 1):
            target = out / f"{run}_{mate}.fastq.gz"
            print(f"{run}\t{mate}\thttps://{url}\t{md5}\t{size}\t{target}")
PY

download_run() {
    local run="$1"
    local row mate url expected_md5 expected_bytes target actual_md5 actual_bytes
    while IFS=$'\t' read -r row mate url expected_md5 expected_bytes target; do
        [[ "$row" == "$run" ]] || continue
        if [[ -s "$target" ]]; then
            actual_md5=$(md5sum "$target" | awk '{print $1}')
        else
            actual_md5=""
        fi
        if [[ "$actual_md5" != "$expected_md5" ]]; then
            printf '[%s] DOWNLOAD %s mate%s -> %s\n' "$(date -Is)" "$run" "$mate" "$target" >> "$DOWNLOAD_LOG"
            wget -c --tries=20 --timeout=60 --retry-connrefused --no-verbose "$url" -O "$target" >> "$DOWNLOAD_LOG" 2>&1
        fi
        actual_bytes=$(stat -c '%s' "$target")
        actual_md5=$(md5sum "$target" | awk '{print $1}')
        if [[ "$actual_bytes" != "$expected_bytes" || "$actual_md5" != "$expected_md5" ]]; then
            printf '[%s] FAIL %s mate%s bytes=%s/%s md5=%s/%s\n' \
                "$(date -Is)" "$run" "$mate" "$actual_bytes" "$expected_bytes" "$actual_md5" "$expected_md5" >> "$DOWNLOAD_LOG"
            return 1
        fi
        printf '[%s] PASS %s mate%s bytes=%s md5=%s\n' \
            "$(date -Is)" "$run" "$mate" "$actual_bytes" "$actual_md5" >> "$DOWNLOAD_LOG"
    done < <(tail -n +2 "$MANIFEST")
}

export -f download_run
export MANIFEST DOWNLOAD_LOG
cut -f1 "$MANIFEST" | tail -n +2 | sort -u | xargs -n1 -P "$JOBS" bash -c 'download_run "$1"' _

find "$OUT" -maxdepth 1 -type f -name '*.fastq.gz' -print0 | xargs -0 -n1 -P 6 gzip -t
"$SEQKIT" stats -T -j 16 "$OUT"/*.fastq.gz > "$AUDIT/rna_fastq_seqkit_stats.tsv"

python3 - "$MANIFEST" "$AUDIT/rna_fastq_seqkit_stats.tsv" > "$STATUS" <<'PY'
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

manifest, stats = Path(sys.argv[1]), Path(sys.argv[2])
stat_rows = {}
with stats.open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        stat_rows[str(Path(row["file"]).resolve())] = row

rows = []
with manifest.open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = Path(row["output"])
        digest = hashlib.md5()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                digest.update(block)
        stats_row = stat_rows.get(str(path.resolve()), {})
        rows.append({
            **row,
            "actual_bytes": path.stat().st_size,
            "actual_md5": digest.hexdigest(),
            "num_seqs": stats_row.get("num_seqs", ""),
        })

by_run = defaultdict(list)
for row in rows:
    by_run[row["run"]].append(row)

print("run\tmate1_num_seqs\tmate2_num_seqs\tpaired_counts_equal\tbytes_md5_pass\tstatus")
for run in sorted(by_run):
    mates = sorted(by_run[run], key=lambda x: int(x["mate"]))
    counts_equal = len(mates) == 2 and mates[0]["num_seqs"] == mates[1]["num_seqs"] and mates[0]["num_seqs"] != ""
    checks = all(
        str(row["actual_bytes"]) == row["expected_bytes"] and row["actual_md5"] == row["expected_md5"]
        for row in mates
    )
    status = "PASS" if counts_equal and checks else "FAIL"
    print(f"{run}\t{mates[0]['num_seqs']}\t{mates[1]['num_seqs']}\t{str(counts_equal).upper()}\t{str(checks).upper()}\t{status}")
    if status != "PASS":
        raise SystemExit(f"Input-integrity gate failed for {run}")
PY

(cd "$OUT" && sha256sum *.fastq.gz | sort > "$AUDIT/rna_fastq_files.sha256")
chmod 0444 "$OUT"/*.fastq.gz
printf '[%s] ALL 9 RUNS PASSED MD5, gzip and paired-read-count gates\n' "$(date -Is)" | tee -a "$DOWNLOAD_LOG"

