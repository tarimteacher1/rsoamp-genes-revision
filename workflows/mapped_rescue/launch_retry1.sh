#!/usr/bin/env bash
set -euo pipefail

root="/path/to/rsoamp/r1_m6_mapped_rescue_20260905"
run_dir="$root/run_retry1"
mkdir -p "$run_dir"
cd "$run_dir"

exec nice -n 10 /path/to/user-home/tools/bpp-workflow-bin/snakemake \
  --snakefile ../workflow/Snakefile \
  --configfile ../input/config.json \
  --cores 8 \
  --printshellcmds \
  --rerun-incomplete \
  --show-failed-logs \
  --latency-wait 60 \
  > master.stdout.log 2>&1
