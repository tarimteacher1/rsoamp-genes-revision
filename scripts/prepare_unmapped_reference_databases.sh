#!/usr/bin/env bash
set -euo pipefail

REVISION=${REVISION:-/path/to/rsoamp/revision_R1_20260902}
REF=${REVISION}/04_unmapped_reads/reference_databases
LOG=${REVISION}/04_unmapped_reads/logs
mkdir -p "${REF}" "${LOG}"

SILVA_BASE=https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports
SSU=SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
LSU=SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz

for file in "${SSU}" "${LSU}"; do
  if [[ ! -s "${REF}/${file}" ]]; then
    curl -L --fail --retry 5 \
      "${SILVA_BASE}/${file}" -o "${REF}/${file}.part"
    mv "${REF}/${file}.part" "${REF}/${file}"
  fi
  curl -L --fail --retry 5 \
    "${SILVA_BASE}/${file}.md5" -o "${REF}/${file}.md5"
  (cd "${REF}" && md5sum -c "${file}.md5")
done

if [[ ! -s "${REF}/SILVA_138.2_SSU_LSU_NR99.fasta" ]]; then
  zcat "${REF}/${SSU}" "${REF}/${LSU}" > "${REF}/SILVA_138.2_SSU_LSU_NR99.fasta.tmp"
  mv "${REF}/SILVA_138.2_SSU_LSU_NR99.fasta.tmp" "${REF}/SILVA_138.2_SSU_LSU_NR99.fasta"
fi

fetch_ncbi_fasta() {
  local accession=$1 output=$2
  if [[ ! -s "${output}" ]]; then
    curl -L --fail --retry 5 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text" \
      -o "${output}.part"
    grep -q '^>' "${output}.part"
    mv "${output}.part" "${output}"
  fi
}

fetch_ncbi_fasta NC_041273 "${REF}/Reaumuria_songarica_chloroplast_NC_041273.fasta"
fetch_ncbi_fasta MW971331 "${REF}/Myricaria_laxiflora_mitochondrion_MW971331_proxy.fasta"

cat > "${REF}/reference_manifest.tsv" <<EOF
reference	accession_or_release	scope	classification_role	limitation
SILVA SSU Ref NR99	138.2	rRNA across cellular life	rRNA-derived read screen	K-mer screen may miss highly divergent or short rRNA fragments
SILVA LSU Ref NR99	138.2	rRNA across cellular life	rRNA-derived read screen	K-mer screen may miss highly divergent or short rRNA fragments
Reaumuria songarica chloroplast	NC_041273	same species	plastid-derived read screen	None beyond reference divergence among individuals
Myricaria laxiflora mitochondrion	MW971331	same family Tamaricaceae	mitochondrial-like read proxy screen	Not a species-specific R. soongarica mitochondrial reference; counts are conservative proxy matches
EOF

sha256sum \
  "${REF}/${SSU}" "${REF}/${LSU}" \
  "${REF}/SILVA_138.2_SSU_LSU_NR99.fasta" \
  "${REF}/Reaumuria_songarica_chloroplast_NC_041273.fasta" \
  "${REF}/Myricaria_laxiflora_mitochondrion_MW971331_proxy.fasta" \
  > "${REF}/reference_databases.sha256"
date --iso-8601=seconds > "${REF}/REFERENCE_DATABASES_COMPLETE.PASS"
echo "PASS unmapped-read reference databases"
