#!/usr/bin/env bash
set -euo pipefail

REVISION=${REVISION:-/path/to/rsoamp/revision_R1_20260902}
MAFFT=${MAFFT:-/path/to/user-home/tools/miniconda3/envs/busco1/bin/mafft}
TRIMAL=${TRIMAL:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/trimal}
IQTREE=${IQTREE:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/iqtree3}
THREADS=${THREADS:-16}
SEED=${SEED:-20260902}
MODEL_SET=${MODEL_SET:-LG,WAG,JTT,VT,DAYHOFF}
RATE_SET=${RATE_SET:-E,I,G4,R4}
PHYLO=${REVISION}/02_phylogeny
mkdir -p "${PHYLO}/logs" "${REVISION}/08_reproducibility/commands"

{ "${MAFFT}" --version 2>&1 || true; } | head -1 > "${PHYLO}/software_versions.txt"
{ "${TRIMAL}" --version 2>&1 || true; } | head -1 >> "${PHYLO}/software_versions.txt"
{ "${IQTREE}" --version 2>&1 || true; } | head -2 >> "${PHYLO}/software_versions.txt"

for family in Defensin Snakin_GASA nsLTP; do
  family_dir=${PHYLO}/${family}
  for mode in full core; do
    input=${family_dir}/${family}.${mode}.faa
    aligned=${family_dir}/${family}.${mode}.mafft_linsi.faa
    trimmed=${family_dir}/${family}.${mode}.mafft_linsi.gappyout.faa
    [[ -s "${input}" ]] || { echo "Missing ${input}" >&2; exit 1; }

    "${MAFFT}" --localpair --maxiterate 1000 --thread "${THREADS}" "${input}" \
      > "${aligned}" 2> "${PHYLO}/logs/${family}.${mode}.mafft.log"
    "${TRIMAL}" -in "${aligned}" -out "${trimmed}" -gappyout \
      > "${PHYLO}/logs/${family}.${mode}.trimal.stdout.log" \
      2> "${PHYLO}/logs/${family}.${mode}.trimal.stderr.log"

    for variant in untrimmed gappyout; do
      if [[ "${variant}" == untrimmed ]]; then alignment=${aligned}; else alignment=${trimmed}; fi
      prefix=${family_dir}/${family}.${mode}.${variant}
      cat > "${REVISION}/08_reproducibility/commands/${family}.${mode}.${variant}.iqtree.command.txt" <<EOF
${IQTREE} -s ${alignment} -m MFP -mset ${MODEL_SET} -mrate ${RATE_SET} -alrt 1000 -bb 1000 -bnni -seed ${SEED} -nt AUTO -ntmax ${THREADS} -pre ${prefix} -redo
EOF
      "${IQTREE}" -s "${alignment}" -m MFP -mset "${MODEL_SET}" -mrate "${RATE_SET}" \
        -alrt 1000 -bb 1000 -bnni \
        -seed "${SEED}" -nt AUTO -ntmax "${THREADS}" -pre "${prefix}" -redo \
        > "${PHYLO}/logs/${family}.${mode}.${variant}.iqtree.stdout.log" \
        2> "${PHYLO}/logs/${family}.${mode}.${variant}.iqtree.stderr.log"
      [[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]] || {
        echo "Incomplete IQ-TREE output: ${prefix}" >&2
        exit 1
      }
    done
  done
done

find "${PHYLO}" -type f \( -name '*.treefile' -o -name '*.contree' -o -name '*.faa' \) \
  -print0 | sort -z | xargs -0 sha256sum > "${PHYLO}/phylogeny_outputs.sha256"
date --iso-8601=seconds > "${PHYLO}/PHYLOGENY_SENSITIVITY_COMPLETE.PASS"
echo "PASS phylogeny sensitivity analyses"
