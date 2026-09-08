#!/usr/bin/env bash
set -u

ROOT=${ROOT:-/path/to/rsoamp/r1_m6_full_read_validation_20260905}
cd "${ROOT}"
nice -n 10 bash workflow/run_isoseq_validation.sh > isoseq.launch.stdout.log 2>&1
code=$?
printf '%s\n' "${code}" > isoseq.run.exitcode
exit "${code}"
