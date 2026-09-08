#!/usr/bin/env bash
set +e
export RUN_ROOT=${RUN_ROOT:-$(pwd)}
nice -n 10 bash workflow/run_overlap_family_audit.sh > launch.stdout.log 2>&1
rc=$?
printf '%s\n' "${rc}" > run.exitcode
exit "${rc}"
