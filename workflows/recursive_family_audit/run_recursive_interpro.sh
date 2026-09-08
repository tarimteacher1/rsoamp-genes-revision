#!/usr/bin/env bash
set -euo pipefail

PROJECT=${PROJECT:-/path/to/rsoamp}
RUN_ROOT=${RUN_ROOT:-${PROJECT}/r1_m6_recursive_family_audit_20260905}
PREP=${PREP:-${RUN_ROOT}/prep}
OUT=${OUT:-${RUN_ROOT}/interpro}
INTERPROSCAN=${INTERPROSCAN:-/path/to/data/tools/software/interproscan-5.77-108.0/interproscan.sh}
JAVA_BIN=${JAVA_BIN:-/path/to/user-home/tools/miniconda3/lib/jvm/bin}
QUERY=${PREP}/recursive_pf00234_start_stop.faa

mkdir -p "${OUT}/logs"
exec > >(tee -a "${OUT}/logs/interproscan.master.log") 2>&1
[[ -s "${QUERY}" ]] || { echo "Missing query: ${QUERY}" >&2; exit 1; }
[[ -x "${INTERPROSCAN}" ]] || { echo "Missing InterProScan: ${INTERPROSCAN}" >&2; exit 1; }
export PATH="${JAVA_BIN}:/usr/bin:/bin"

"${INTERPROSCAN}" -i "${QUERY}" -b "${OUT}/recursive_pf00234.interproscan" \
  -f TSV,GFF3 -cpu 1 -dp -iprlookup -goterms
[[ -s "${OUT}/recursive_pf00234.interproscan.tsv" ]] || { echo "Missing InterProScan TSV" >&2; exit 1; }
[[ -s "${OUT}/recursive_pf00234.interproscan.gff3" ]] || { echo "Missing InterProScan GFF3" >&2; exit 1; }
sha256sum "${QUERY}" "${OUT}/recursive_pf00234.interproscan.tsv" \
  "${OUT}/recursive_pf00234.interproscan.gff3" > "${OUT}/interproscan.sha256"
date --iso-8601=seconds > "${OUT}/RECURSIVE_INTERPRO_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS recursive InterProScan"
