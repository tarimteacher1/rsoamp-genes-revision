#!/usr/bin/env bash
set -u

ROOT=${ROOT:-/path/to/rsoamp/r1_m6_recursive_family_audit_20260905}
cd "${ROOT}"
nice -n 10 bash workflow/run_recursive_family_audit.sh > launch.stdout.log 2>&1
code=$?
printf '%s\n' "${code}" > run.exitcode
exit "${code}"
