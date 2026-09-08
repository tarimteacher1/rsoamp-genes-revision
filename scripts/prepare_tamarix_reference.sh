#!/usr/bin/env bash
# Acquire and screen the closest chromosome-scale Tamaricaceae reference used in revision R1.
set -euo pipefail

REV=/path/to/rsoamp/revision_R1_20260902
ORIG=/path/to/rsoamp
OUT=$REV/05_comparative_genomics
IN=$OUT/inputs/Tamarix_austromongolica
LOG=$OUT/logs
PRED=$OUT/predictions
mkdir -p "$IN" "$LOG" "$PRED"

PEP=$IN/tau.longest_pep.fasta
GFF=$IN/tau.longest.gff3

download_verify() {
    local url=$1 expected=$2 dest=$3
    if [[ -s "$dest" ]]; then
        local existing
        existing=$(md5sum "$dest" | cut -d' ' -f1)
        if [[ "$existing" == "$expected" ]]; then
            echo "Reusing verified file: $dest"
            return 0
        fi
        chmod u+w "$dest"
    fi
    wget -c --tries=20 --timeout=60 --retry-connrefused --no-verbose "$url" -O "$dest"
    local observed
    observed=$(md5sum "$dest" | cut -d' ' -f1)
    if [[ "$observed" != "$expected" ]]; then
        echo "MD5 mismatch for $dest: expected=$expected observed=$observed" >&2
        exit 1
    fi
}

download_verify https://ndownloader.figshare.com/files/44299196 7121178880871c51cee457be252fb48e "$PEP"
download_verify https://ndownloader.figshare.com/files/44299208 fb6450acc0e9cd8ca6b938aaf2a5ce9d "$GFF"

cat > "$OUT/reference_source_manifest.tsv" <<'EOF'
species	taxonomy	assembly_accession	assembly_level	publication_doi	data_doi	file	url	expected_md5	license
Tamarix austromongolica	Tamaricaceae; Caryophyllales	GCA_039764185.1	Chromosome	10.1093/dnares/dsae021	10.6084/m9.figshare.25106726.v1	tau.longest_pep.fasta	https://ndownloader.figshare.com/files/44299196	7121178880871c51cee457be252fb48e	CC BY 4.0
Tamarix austromongolica	Tamaricaceae; Caryophyllales	GCA_039764185.1	Chromosome	10.1093/dnares/dsae021	10.6084/m9.figshare.25106726.v1	tau.longest.gff3	https://ndownloader.figshare.com/files/44299208	fb6450acc0e9cd8ca6b938aaf2a5ce9d	CC BY 4.0
EOF

SEQKIT=/path/to/user-home/tools/miniconda3/envs/cotton/bin/seqkit
HMMSEARCH=/path/to/user-home/tools/miniconda3/envs/busco1/bin/hmmsearch
"$SEQKIT" stats -T "$PEP" > "$OUT/tamarix_proteome_stats.tsv"

declare -A MODELS=(
    [Defensin]=PF00304
    [Thionin]=PF00321
    [Snakin_GASA]=PF02704
    [nsLTP]=PF00234
    [Hevein_like]=PF00187
    [Cyclotide]=PF03784
)

for family in Defensin Thionin Snakin_GASA nsLTP Hevein_like Cyclotide; do
    pfam=${MODELS[$family]}
    "$HMMSEARCH" --noali -E 1e-5 --domE 1e-5 \
        --tblout "$PRED/${family}.tblout" \
        --domtblout "$PRED/${family}.domtblout" \
        "$ORIG/data/hmm/${pfam}.hmm" "$PEP" \
        > "$LOG/${family}.hmmsearch.log"
done

python3 - "$PRED" "$PEP" <<'PY'
import csv
import sys
from pathlib import Path

pred = Path(sys.argv[1])
pep = Path(sys.argv[2])
families = ["Defensin", "Thionin", "Snakin_GASA", "nsLTP", "Hevein_like", "Cyclotide"]
hits = {}
rows = []
for family in families:
    path = pred / f"{family}.tblout"
    ids = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split()
            target = fields[0]
            evalue = fields[4]
            score = fields[5]
            ids.add(target)
            rows.append({"protein_id": target, "family_hmm": family, "full_evalue": evalue, "full_score": score})
    hits[family] = ids

with (pred / "tamarix_amp_hmm_hits.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["protein_id", "family_hmm", "full_evalue", "full_score"], delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda row: (row["family_hmm"], row["protein_id"])))

wanted = set().union(*hits.values())
seqs = {}
current = None
for raw in pep.open():
    line = raw.strip()
    if line.startswith(">"):
        current = line[1:].split()[0]
        if current in wanted:
            seqs[current] = []
    elif current in seqs:
        seqs[current].append(line)
with (pred / "tamarix_amp_hmm_candidates.faa").open("w") as handle:
    for protein_id in sorted(seqs):
        handle.write(f">{protein_id}\n{''.join(seqs[protein_id])}\n")

with (pred / "tamarix_amp_hmm_counts.tsv").open("w") as handle:
    handle.write("family_hmm\tn_hits\n")
    for family in families:
        handle.write(f"{family}\t{len(hits[family])}\n")
    handle.write(f"union\t{len(wanted)}\n")
print("Tamarix HMM candidate union:", len(wanted))
PY

sha256sum "$PEP" "$GFF" "$PRED"/*.tblout "$PRED"/*.domtblout "$PRED/tamarix_amp_hmm_candidates.faa" \
    > "$OUT/input_and_hmm_outputs.sha256"

chmod a-w "$PEP" "$GFF"
echo "Tamarix reference preparation PASS"
