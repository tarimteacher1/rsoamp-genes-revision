#!/usr/bin/env bash
set -u

ROOT=${ROOT:-/path/to/rsoamp/r1_m6_full_read_validation_20260905}
cd "${ROOT}"
nice -n 10 bash workflow/run_full_read_validation.sh > launch.stdout.log 2>&1
code=$?
printf '%s\n' "${code}" > run.exitcode
exit "${code}"
