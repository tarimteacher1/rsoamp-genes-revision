#!/usr/bin/env bash
set -euo pipefail

PROJECT=${PROJECT:-/path/to/rsoamp}
REVISION=${REVISION:-${PROJECT}/revision_R1_20260902}
SOURCE_RESULTS=${SOURCE_RESULTS:-${PROJECT}/r1_m6_mapped_rescue_20260905/run_retry1/results}
RUN_ROOT=${RUN_ROOT:-${PROJECT}/r1_m6_full_read_validation_20260905}
PREP=${PREP:-${RUN_ROOT}/prep}
OUT=${OUT:-${RUN_ROOT}/isoseq}
THREADS=${THREADS:-4}

MINIMAP2=${MINIMAP2:-${REVISION}/08_reproducibility/envs/isoseq_mapping/bin/minimap2}
SAMTOOLS=${SAMTOOLS:-/path/to/user-home/tools/miniconda3/envs/as_splice/bin/samtools}
SEQKIT=${SEQKIT:-/path/to/user-home/tools/miniconda3/envs/seqkit/bin/seqkit}
STRINGTIE=${STRINGTIE:-/path/to/user-home/tools/miniconda3/envs/genome_annot/bin/stringtie}
GFFCOMPARE=${GFFCOMPARE:-/path/to/user-home/tools/miniconda3/envs/genome_annot/bin/gffcompare}
GFFREAD=${GFFREAD:-/path/to/user-home/tools/miniconda3/envs/commom_tools2/bin/gffread}
GETORF=${GETORF:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/getorf}
DIAMOND=${DIAMOND:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/diamond}
HMMSEARCH=${HMMSEARCH:-/path/to/user-home/tools/miniconda3/envs/cotton/bin/hmmsearch}

INPUT=${REVISION}/03_annotation_rescue/isoseq/SRR27540882.sra_derived.fastq.gz
GENOME=${PROJECT}/data/genome.fa
REFERENCE_GTF=${REVISION}/04_unmapped_reads/reference_annotation.gtf
WINDOWS=${PREP}/candidate_windows_10kb.bed
TARGET_FASTA=${OUT}/candidate_windows_10kb.fasta
SCREEN_BAM=${OUT}/local_window_screen.bam
NAMES=${OUT}/candidate_read_names.txt
CANDIDATE_FASTQ=${OUT}/candidate_reads.fastq.gz
GLOBAL_BAM=${OUT}/candidate_reads.global_remap.bam

mkdir -p "${OUT}/logs" "${OUT}/screen"
exec > >(tee -a "${OUT}/logs/isoseq_validation.master.log") 2>&1
for path in "${INPUT}" "${GENOME}" "${REFERENCE_GTF}" "${WINDOWS}" "${PREP}/intergenic_candidate_models.gtf"; do
  [[ -s "${path}" ]] || { echo "Missing required input: ${path}" >&2; exit 1; }
done
for tool in "${MINIMAP2}" "${SAMTOOLS}" "${SEQKIT}" "${STRINGTIE}" "${GFFCOMPARE}" "${GFFREAD}" "${GETORF}" "${DIAMOND}" "${HMMSEARCH}"; do
  [[ -x "${tool}" ]] || { echo "Missing executable: ${tool}" >&2; exit 1; }
done

if [[ ! -s "${TARGET_FASTA}" ]]; then
  : > "${TARGET_FASTA}"
  while IFS=$'\t' read -r chrom start0 end0 name rest; do
    "${SAMTOOLS}" faidx "${GENOME}" "${chrom}:$((start0 + 1))-${end0}" >> "${TARGET_FASTA}"
  done < "${WINDOWS}"
fi

if [[ ! -s "${SCREEN_BAM}" ]]; then
  echo "[$(date --iso-8601=seconds)] START Iso-Seq local-window screen"
  set +e
  "${MINIMAP2}" -ax splice:hq -uf --secondary=yes -N 20 --MD -t "${THREADS}" \
    "${TARGET_FASTA}" "${INPUT}" 2> "${OUT}/logs/local_window_screen.minimap2.log" \
    | "${SAMTOOLS}" view -bh -F 4 -o "${SCREEN_BAM}.tmp" -
  statuses=("${PIPESTATUS[@]}")
  set -e
  if [[ "${statuses[0]}" -ne 0 || "${statuses[1]}" -ne 0 ]]; then
    echo "Iso-Seq local screen failed: minimap2=${statuses[0]} samtools=${statuses[1]}" >&2
    exit 2
  fi
  mv "${SCREEN_BAM}.tmp" "${SCREEN_BAM}"
fi
"${SAMTOOLS}" quickcheck -v "${SCREEN_BAM}"
"${SAMTOOLS}" view "${SCREEN_BAM}" | cut -f1 | LC_ALL=C sort -u > "${NAMES}"
if [[ -s "${NAMES}" && ! -s "${CANDIDATE_FASTQ}" ]]; then
  "${SEQKIT}" grep -f "${NAMES}" "${INPUT}" -o "${CANDIDATE_FASTQ}.tmp" \
    > "${OUT}/logs/seqkit_extract.stdout.log" 2> "${OUT}/logs/seqkit_extract.stderr.log"
  mv "${CANDIDATE_FASTQ}.tmp" "${CANDIDATE_FASTQ}"
elif [[ ! -s "${NAMES}" ]]; then
  gzip -c </dev/null > "${CANDIDATE_FASTQ}"
fi

if [[ -s "${NAMES}" && ! -s "${GLOBAL_BAM}" ]]; then
  echo "[$(date --iso-8601=seconds)] START Iso-Seq whole-genome remap of candidate reads"
  set +e
  "${MINIMAP2}" -ax splice:hq -uf --secondary=yes -N 20 --MD -t "${THREADS}" \
    "${GENOME}" "${CANDIDATE_FASTQ}" 2> "${OUT}/logs/global_candidate_remap.minimap2.log" \
    | "${SAMTOOLS}" view -bh -L "${WINDOWS}" - \
    | "${SAMTOOLS}" sort -@ 2 -m 1G -o "${GLOBAL_BAM}.tmp" -
  statuses=("${PIPESTATUS[@]}")
  set -e
  if [[ "${statuses[0]}" -ne 0 || "${statuses[1]}" -ne 0 || "${statuses[2]}" -ne 0 ]]; then
    echo "Iso-Seq global remap failed: ${statuses[*]}" >&2
    exit 2
  fi
  mv "${GLOBAL_BAM}.tmp" "${GLOBAL_BAM}"
  "${SAMTOOLS}" index -@ 2 "${GLOBAL_BAM}"
fi

if [[ -s "${GLOBAL_BAM}" ]]; then
  "${SAMTOOLS}" quickcheck -v "${GLOBAL_BAM}"
  [[ -s "${GLOBAL_BAM}.bai" ]] || "${SAMTOOLS}" index -@ 2 "${GLOBAL_BAM}"
  "${SAMTOOLS}" flagstat -O json "${GLOBAL_BAM}" > "${OUT}/global_remap.flagstat.json"
  "${SAMTOOLS}" bedcov "${WINDOWS}" "${GLOBAL_BAM}" > "${OUT}/candidate_windows.bedcov.tsv"
  "${STRINGTIE}" "${GLOBAL_BAM}" -L -p "${THREADS}" -c 1 -s 1 -G "${REFERENCE_GTF}" \
    -o "${OUT}/screen/isoseq_candidate_models.gtf" 2> "${OUT}/logs/stringtie_longread.log"
  "${GFFCOMPARE}" -r "${PREP}/intergenic_candidate_models.gtf" -o "${OUT}/screen/vs_subset_candidate" \
    "${OUT}/screen/isoseq_candidate_models.gtf"
  "${GFFCOMPARE}" -r "${REFERENCE_GTF}" -o "${OUT}/screen/vs_reference" \
    "${OUT}/screen/isoseq_candidate_models.gtf"
  "${GFFREAD}" "${OUT}/screen/isoseq_candidate_models.gtf" -g "${GENOME}" \
    -w "${OUT}/screen/isoseq_candidate_transcripts.fna"
  PATH="$(dirname "${GETORF}"):${PATH}" "${GETORF}" \
    -sequence "${OUT}/screen/isoseq_candidate_transcripts.fna" \
    -outseq "${OUT}/screen/isoseq_candidate_orfs.faa" -find 0 -minsize 90 -reverse N -table 0 -auto
  "${DIAMOND}" blastp --query "${OUT}/screen/isoseq_candidate_orfs.faa" \
    --db "${SOURCE_RESULTS}/screen/references.dmnd" --out "${OUT}/screen/isoseq_orf_reference_hits.tsv" \
    --sensitive --evalue 0.001 --id 20 --max-target-seqs 0 --threads "${THREADS}" \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp scovhsp evalue bitscore
  "${HMMSEARCH}" --cpu "${THREADS}" --cut_ga --noali \
    --domtblout "${OUT}/screen/isoseq_orf_family_domains.domtbl" \
    "${PROJECT}/data/hmm/AMP_models.hmm" "${OUT}/screen/isoseq_candidate_orfs.faa" \
    > "${OUT}/logs/isoseq_hmmsearch.log"
fi

{
  printf 'metric\tvalue\tinterpretation\n'
  printf 'total_isoseq_reads\t%s\tVerified public Iso-Seq FASTQ records\n' \
    "$("${SEQKIT}" stats -T "${INPUT}" | awk 'NR==2 {print $4}')"
  printf 'local_window_candidate_reads\t%s\tReads with at least one local-window alignment before whole-genome remapping\n' \
    "$(wc -l < "${NAMES}")"
  printf 'global_remap_reference\twhole_1.28_Gb_genome\tCandidate reads were remapped globally before target-window filtering\n'
} > "${OUT}/isoseq_mapping_metrics.tsv"

sha256sum "${INPUT}" "${WINDOWS}" "${TARGET_FASTA}" "${NAMES}" "${CANDIDATE_FASTQ}" \
  > "${OUT}/input_output_identities.sha256"
date --iso-8601=seconds > "${OUT}/ISOSEQ_CANDIDATE_VALIDATION_RAW_COMPLETE.PASS"
echo "[$(date --iso-8601=seconds)] PASS Iso-Seq candidate validation raw workflow"
