#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/rsoamp/revision_R1_20260902}"
FASTP="${FASTP:-/path/to/user-home/tools/miniconda3/envs/salmon_env/bin/fastp}"
THREADS_PER_JOB="${THREADS_PER_JOB:-12}"
MAX_JOBS="${MAX_JOBS:-2}"
RAW="$BASE/04_unmapped_reads/raw_fastq"
OUT="$BASE/04_unmapped_reads/fastp_qc"
SAMPLES="$BASE/04_unmapped_reads/rnaseq_sample_sheet.tsv"
mkdir -p "$OUT"

[[ -x "$FASTP" ]] || { echo "Missing fastp executable: $FASTP" >&2; exit 2; }
[[ -s "$SAMPLES" ]] || { echo "Missing sample sheet: $SAMPLES" >&2; exit 2; }

qc_one() {
    local run=$1
    local r1="$RAW/${run}_1.fastq.gz"
    local r2="$RAW/${run}_2.fastq.gz"
    local json="$OUT/${run}.fastp.json"
    local html="$OUT/${run}.fastp.html"
    local log="$OUT/${run}.fastp.stderr.log"
    if [[ -s "$OUT/${run}.FASTP_QC_COMPLETE.PASS" ]]; then
        return 0
    fi
    "$FASTP" -i "$r1" -I "$r2" --stdout --detect_adapter_for_pe \
        --thread "$THREADS_PER_JOB" --json "$json" --html "$html" \
        > /dev/null 2> "$log"
    [[ -s "$json" && -s "$html" ]] || { echo "fastp QC failed for $run" >&2; return 1; }
    date --iso-8601=seconds > "$OUT/${run}.FASTP_QC_COMPLETE.PASS"
}
export -f qc_one
export RAW OUT FASTP THREADS_PER_JOB
while IFS=$'\t' read -r run treatment replicate; do
    qc_one "$run" &
    while (( $(jobs -pr | wc -l) >= MAX_JOBS )); do wait -n; done
done < <(tail -n +2 "$SAMPLES")
wait

python3 - "$OUT" "$SAMPLES" <<'PY'
import csv, json, sys
from pathlib import Path
out, samples = map(Path, sys.argv[1:])
rows=[]
with samples.open() as handle:
    for sample in csv.DictReader(handle, delimiter='\t'):
        report=json.loads((out/f"{sample['run']}.fastp.json").read_text())
        before=report['summary']['before_filtering']; after=report['summary']['after_filtering']
        filtering=report.get('filtering_result',{}); adapter=report.get('adapter_cutting',{})
        rows.append({
            **sample,
            'before_total_reads':before.get('total_reads',''),
            'before_total_bases':before.get('total_bases',''),
            'before_q20_rate':before.get('q20_rate',''),
            'before_q30_rate':before.get('q30_rate',''),
            'before_gc_content':before.get('gc_content',''),
            'after_total_reads':after.get('total_reads',''),
            'after_q30_rate':after.get('q30_rate',''),
            'passed_filter_reads':filtering.get('passed_filter_reads',''),
            'low_quality_reads':filtering.get('low_quality_reads',''),
            'too_many_N_reads':filtering.get('too_many_N_reads',''),
            'too_short_reads':filtering.get('too_short_reads',''),
            'adapter_trimmed_reads':adapter.get('adapter_trimmed_reads',''),
            'qc_interpretation':'QC simulation only; original verified FASTQs were retained and used for quantification',
        })
with (out.parent/'per_sample_fastp_qc.tsv').open('w',newline='') as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
print('PASS fastp QC',len(rows),'samples')
PY

[[ $(find "$OUT" -name '*.FASTP_QC_COMPLETE.PASS' | wc -l) -eq 9 ]] || exit 4
date --iso-8601=seconds > "$BASE/04_unmapped_reads/FASTP_QC_COMPLETE.PASS"
