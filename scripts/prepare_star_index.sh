#!/usr/bin/env bash
set -euo pipefail

ORIGINAL=${ORIGINAL:-/path/to/rsoamp}
REVISION=${REVISION:-${ORIGINAL}/revision_R1_20260902}
STAR=${STAR:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/STAR}
GFFREAD=${GFFREAD:-/path/to/user-home/tools/miniconda3/envs/commom_tools2/bin/gffread}
THREADS=${THREADS:-48}

WORK=${REVISION}/04_unmapped_reads
INDEX=${WORK}/star_index_sjdb149
GTF=${WORK}/reference_annotation.gtf
LOGDIR=${WORK}/logs
mkdir -p "${INDEX}" "${LOGDIR}" "${REVISION}/08_reproducibility/commands"

exec > >(tee "${LOGDIR}/prepare_star_index.log") 2>&1
echo "[$(date --iso-8601=seconds)] START STAR index preparation"
"${STAR}" --version
"${GFFREAD}" --version 2>&1 | head -1 || true

if [[ ! -s "${GTF}" ]]; then
  "${GFFREAD}" "${ORIGINAL}/data/genome.gff" -T -o "${GTF}"
fi

cat > "${REVISION}/08_reproducibility/commands/star_genome_generate.command.txt" <<EOF
${STAR} --runThreadN ${THREADS} --runMode genomeGenerate --genomeDir ${INDEX} --genomeFastaFiles ${ORIGINAL}/data/genome.fa --sjdbGTFfile ${GTF} --sjdbOverhang 149 --genomeSAindexNbases 14
EOF

"${STAR}" \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${INDEX}" \
  --genomeFastaFiles "${ORIGINAL}/data/genome.fa" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 14

required=(Genome SA SAindex chrLength.txt chrName.txt chrNameLength.txt genomeParameters.txt sjdbInfo.txt)
for filename in "${required[@]}"; do
  [[ -s "${INDEX}/${filename}" ]] || {
    echo "FAIL missing or empty STAR index component: ${filename}" >&2
    exit 1
  }
done

{
  printf 'file\tsize_bytes\n'
  for filename in "${required[@]}"; do
    printf '%s\t%s\n' "${filename}" "$(stat -c %s "${INDEX}/${filename}")"
  done
} > "${WORK}/star_index_integrity.tsv"

date --iso-8601=seconds > "${INDEX}/INDEX_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS STAR index preparation"
