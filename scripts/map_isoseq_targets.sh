#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/rsoamp/revision_R1_20260902}"
PROJECT="${PROJECT:-/path/to/rsoamp}"
THREADS="${THREADS:-32}"
MINIMAP2="${MINIMAP2:-$BASE/08_reproducibility/envs/isoseq_mapping/bin/minimap2}"
SAMTOOLS="${SAMTOOLS:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/samtools}"
SEQKIT="${SEQKIT:-/path/to/user-home/tools/miniconda3/envs/seqkit/bin/seqkit}"

INPUT="$BASE/03_annotation_rescue/isoseq/SRR27540882.sra_derived.fastq.gz"
GENOME="$PROJECT/data/genome.fa"
LOCI="$BASE/03_annotation_rescue/backfill_locus_audit_preliminary.tsv"
OUTDIR="$BASE/03_annotation_rescue/isoseq_target_alignment"
LOGDIR="$BASE/03_annotation_rescue/logs"
WINDOWS="$OUTDIR/backfill_windows_10kb.bed"
CODING="$OUTDIR/backfill_coding_intervals.bed"
BAM="$OUTDIR/SRR27540882.backfill_targets.bam"
TARGET_FASTA="$OUTDIR/backfill_windows_10kb.fasta"
SCREEN_BAM="$OUTDIR/SRR27540882.local_window_screen.bam"
SCREEN_NAMES="$OUTDIR/SRR27540882.local_window_candidate_read_names.txt"
CANDIDATE_FASTQ="$OUTDIR/SRR27540882.local_window_candidates.fastq.gz"
SUMMARY="$BASE/03_annotation_rescue/isoseq_backfill_support.tsv"

mkdir -p "$OUTDIR" "$LOGDIR" "$BASE/08_reproducibility/versions"

for path in "$INPUT" "$GENOME" "$LOCI"; do
    if [[ ! -s "$path" ]]; then
        echo "Missing or empty required input: $path" >&2
        exit 2
    fi
done

if [[ ! -x "$MINIMAP2" || ! -x "$SAMTOOLS" || ! -x "$SEQKIT" ]]; then
    echo "Missing executable: MINIMAP2=$MINIMAP2 SAMTOOLS=$SAMTOOLS SEQKIT=$SEQKIT" >&2
    exit 2
fi

awk -F '\t' 'BEGIN{OFS="\t"}
NR==1 {
  for(i=1;i<=NF;i++) h[$i]=i
  next
}
{
  split($(h["historical_locus"]), a, /[:()\-]/)
  chr=a[1]
  cs=$(h["coding_start"])+0
  ce=$(h["coding_end"])+0
  lo=(cs<ce?cs:ce)
  hi=(cs>ce?cs:ce)
  start=lo-10001
  if(start<0) start=0
  end=hi+10000
  strand=($(h["historical_locus"]) ~ /\(\+\)$/ ? "+" : "-")
  print chr,start,end,$(h["member_id"]),0,strand
}' "$LOCI" | sort -k1,1 -k2,2n > "$WINDOWS"

awk -F '\t' 'BEGIN{OFS="\t"}
NR==1 {
  for(i=1;i<=NF;i++) h[$i]=i
  next
}
{
  split($(h["historical_locus"]), a, /[:()\-]/)
  chr=a[1]
  cs=$(h["coding_start"])+0
  ce=$(h["coding_end"])+0
  lo=(cs<ce?cs:ce)
  hi=(cs>ce?cs:ce)
  strand=($(h["historical_locus"]) ~ /\(\+\)$/ ? "+" : "-")
  print chr,lo-1,hi,$(h["member_id"]),0,strand
}' "$LOCI" | sort -k1,1 -k2,2n > "$CODING"

"$MINIMAP2" --version > "$BASE/08_reproducibility/versions/minimap2.txt"
"$SAMTOOLS" --version | head -n 1 > "$BASE/08_reproducibility/versions/samtools_isoseq.txt"
"$SEQKIT" version > "$BASE/08_reproducibility/versions/seqkit_isoseq.txt" 2>&1

if [[ ! -s "$TARGET_FASTA" ]]; then
    : > "$TARGET_FASTA"
    while IFS=$'\t' read -r chrom start0 end0 member score strand; do
        "$SAMTOOLS" faidx "$GENOME" "${chrom}:$((start0 + 1))-${end0}" >> "$TARGET_FASTA"
    done < "$WINDOWS"
fi

if [[ ! -s "$SCREEN_BAM" ]]; then
    rm -f "$SCREEN_BAM"
    set +e
    "$MINIMAP2" -ax splice:hq -uf --secondary=yes -N 20 --MD -t "$THREADS" \
        "$TARGET_FASTA" "$INPUT" \
        2> "$LOGDIR/SRR27540882.local_window_screen.minimap2.log" \
      | "$SAMTOOLS" view -bh -F 4 -o "$SCREEN_BAM" -
    status=("${PIPESTATUS[@]}")
    set -e
    if [[ ${status[0]} -ne 0 || ${status[1]} -ne 0 ]]; then
        echo "Iso-Seq local-window screening failed: ${status[*]}" >&2
        exit 3
    fi
fi

"$SAMTOOLS" quickcheck -v "$SCREEN_BAM"
"$SAMTOOLS" view "$SCREEN_BAM" | cut -f1 | LC_ALL=C sort -u > "$SCREEN_NAMES"
if [[ -s "$SCREEN_NAMES" ]]; then
    rm -f "$CANDIDATE_FASTQ"
    "$SEQKIT" grep -f "$SCREEN_NAMES" "$INPUT" -o "$CANDIDATE_FASTQ" \
        > "$LOGDIR/SRR27540882.seqkit_candidate_extract.stdout.log" \
        2> "$LOGDIR/SRR27540882.seqkit_candidate_extract.stderr.log"
else
    gzip -c </dev/null > "$CANDIDATE_FASTQ"
fi

if [[ ! -s "$BAM" || ! -s "$BAM.bai" ]]; then
    rm -f "$BAM" "$BAM.bai" "$OUTDIR/SRR27540882.backfill_targets.unsorted.bam"
    set +e
    "$MINIMAP2" -ax splice:hq -uf --secondary=yes -N 5 --MD -t "$THREADS" \
        "$GENOME" "$CANDIDATE_FASTQ" \
        2> "$LOGDIR/SRR27540882.global_candidate_remap.minimap2.log" \
      | "$SAMTOOLS" view -bh -L "$WINDOWS" - \
      | "$SAMTOOLS" sort -@ 8 -m 1G -o "$BAM" -
    status=("${PIPESTATUS[@]}")
    set -e
    if [[ ${status[0]} -ne 0 || ${status[1]} -ne 0 || ${status[2]} -ne 0 ]]; then
        echo "Iso-Seq target alignment pipeline failed: ${status[*]}" >&2
        exit 3
    fi
    "$SAMTOOLS" index -@ 8 "$BAM"
fi

"$SAMTOOLS" quickcheck -v "$BAM"
python3 "$BASE/03_annotation_rescue/summarize_target_alignment_support.py" \
    --loci "$LOCI" \
    --bam "$BAM" \
    --samtools "$SAMTOOLS" \
    --output "$SUMMARY"

{
    sha256sum "$INPUT" "$GENOME" "$LOCI" "$WINDOWS" "$CODING" "$TARGET_FASTA" "$SCREEN_NAMES" "$CANDIDATE_FASTQ" "$BAM"
} > "$BASE/03_annotation_rescue/isoseq_target_alignment.sha256"

{
    printf 'metric\tvalue\tinterpretation\n'
    printf 'total_isoseq_reads\t%s\tVerified SRR27540882 FASTQ records\n' "$("$SEQKIT" stats -T "$INPUT" | awk 'NR==2{print $4}')"
    printf 'local_window_candidate_reads\t%s\tReads with at least one alignment to a six-locus plus-or-minus-10-kb screening window\n' "$(wc -l < "$SCREEN_NAMES")"
    printf 'global_remap_reference\twhole_1.28_Gb_genome\tCandidate reads were remapped to the complete genome before target-window filtering\n'
} > "$BASE/03_annotation_rescue/isoseq_two_stage_mapping_metrics.tsv"

if [[ $(awk 'END{print NR-1}' "$SUMMARY") -ne 6 ]]; then
    echo "Expected six Iso-Seq locus rows in $SUMMARY" >&2
    exit 4
fi

printf 'PASS\tIso-Seq target alignment and six-locus support summary completed\n' \
    > "$BASE/03_annotation_rescue/ISOSEQ_TARGET_ALIGNMENT_COMPLETE.PASS"
cat "$SUMMARY"
