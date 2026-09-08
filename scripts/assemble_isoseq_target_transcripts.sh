#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE:-/path/to/rsoamp/revision_R1_20260902}"
PROJECT="${PROJECT:-/path/to/rsoamp}"
STRINGTIE="${STRINGTIE:-/path/to/user-home/tools/miniconda3/envs/genome_annot/bin/stringtie}"
GFFREAD="${GFFREAD:-/path/to/user-home/tools/miniconda3/envs/starsolo/bin/gffread}"
THREADS="${THREADS:-16}"

DIR="$BASE/03_annotation_rescue/isoseq_target_alignment"
BAM="$DIR/SRR27540882.backfill_targets.bam"
GTF="$DIR/SRR27540882.backfill_targets.stringtie_longread.gtf"
TRANSCRIPTS="$DIR/SRR27540882.backfill_targets.stringtie_longread.fa"
LOG="$BASE/03_annotation_rescue/logs/SRR27540882.stringtie_longread.log"

[[ -s "$BAM" && -s "$BAM.bai" ]] || { echo "Missing target Iso-Seq BAM or index" >&2; exit 2; }
for tool in "$STRINGTIE" "$GFFREAD"; do
    [[ -x "$tool" ]] || { echo "Missing executable: $tool" >&2; exit 2; }
done

"$STRINGTIE" "$BAM" -L -p "$THREADS" -c 1 -s 1 \
    -G "$PROJECT/data/genome.gff" \
    -o "$GTF" \
    2> "$LOG"
[[ -s "$GTF" ]] || { echo "StringTie produced an empty GTF" >&2; exit 3; }

"$GFFREAD" "$GTF" -g "$PROJECT/data/genome.fa" -w "$TRANSCRIPTS" \
    > "$BASE/03_annotation_rescue/logs/SRR27540882.gffread_target_transcripts.stdout.log" \
    2> "$BASE/03_annotation_rescue/logs/SRR27540882.gffread_target_transcripts.stderr.log"
[[ -s "$TRANSCRIPTS" ]] || { echo "gffread produced no target transcripts" >&2; exit 3; }

python3 "$BASE/03_annotation_rescue/evaluate_isoseq_gene_models.py" \
    --historical-peptides "$PROJECT/intermediate/members_final.faa" \
    --backfill-audit "$BASE/03_annotation_rescue/backfill_locus_audit_preliminary.tsv" \
    --transcripts "$TRANSCRIPTS" \
    --gtf "$GTF" \
    --output "$BASE/03_annotation_rescue/isoseq_backfill_gene_models.tsv"

sha256sum "$BAM" "$GTF" "$TRANSCRIPTS" \
    > "$BASE/03_annotation_rescue/isoseq_target_transcript_models.sha256"
printf 'PASS\tIso-Seq target transcript assembly and complete-ORF evaluation finished\n' \
    > "$BASE/03_annotation_rescue/ISOSEQ_TARGET_TRANSCRIPT_ASSEMBLY_COMPLETE.PASS"
cat "$BASE/03_annotation_rescue/isoseq_backfill_gene_models.tsv"
