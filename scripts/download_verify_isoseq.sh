#!/usr/bin/env bash
set -euo pipefail

REVISION=${REVISION:-/path/to/rsoamp/revision_R1_20260902}
SEQKIT=${SEQKIT:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/seqkit}
FASTERQ=${FASTERQ:-/path/to/user-home/tools/miniconda3/envs/sra_tools/bin/fasterq-dump}
VDB_VALIDATE=${VDB_VALIDATE:-/path/to/user-home/tools/miniconda3/envs/sra_tools/bin/vdb-validate}
OUT=${REVISION}/03_annotation_rescue/isoseq
LOG=${REVISION}/03_annotation_rescue/logs
SRA=${OUT}/SRR27540882.sra
FILE=${OUT}/SRR27540882.sra_derived.fastq.gz
URL=https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR27540882/SRR27540882
EXPECTED_SRA_BYTES=15794093706
EXPECTED_READS=17136488
EXPECTED_BASES=62351585955
mkdir -p "${OUT}" "${LOG}"

exec > >(tee "${LOG}/download_verify_SRR27540882.log") 2>&1
echo "[$(date --iso-8601=seconds)] START SRR27540882 Iso-Seq SRA acquisition"
wget -c --tries=30 --timeout=60 --retry-connrefused --no-verbose "${URL}" -O "${SRA}"

observed_bytes=$(stat -c %s "${SRA}")
[[ "${observed_bytes}" -eq "${EXPECTED_SRA_BYTES}" ]] || {
  echo "SRA byte-size mismatch: expected=${EXPECTED_SRA_BYTES} observed=${observed_bytes}" >&2
  exit 1
}
"${VDB_VALIDATE}" "${SRA}"
rm -f "${FILE}"
"${FASTERQ}" --split-spot -Z --skip-technical -e 32 "${SRA}" | pigz -p 16 > "${FILE}"
gzip -t "${FILE}"
"${SEQKIT}" stats -T -a "${FILE}" > "${OUT}/SRR27540882_seqkit_stats.tsv"
python3 - "${OUT}/SRR27540882_seqkit_stats.tsv" "${EXPECTED_READS}" "${EXPECTED_BASES}" <<'PY'
import csv
import sys
row = next(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
observed_reads = int(row["num_seqs"].replace(",", ""))
observed_bases = int(row["sum_len"].replace(",", ""))
expected_reads, expected_bases = map(int, sys.argv[2:])
if (observed_reads, observed_bases) != (expected_reads, expected_bases):
    raise SystemExit(
        f"SRA-derived FASTQ count mismatch: observed={observed_reads},{observed_bases} "
        f"expected={expected_reads},{expected_bases}"
    )
PY
sha256sum "${SRA}" "${FILE}" > "${OUT}/SRR27540882.sha256"

date --iso-8601=seconds > "${OUT}/DOWNLOAD_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS SRR27540882 SRA validation, gzip and ENA run-count gates"
