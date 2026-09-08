"""Evidence-bounded supplement for genome-mapped, transcript-unassigned reads."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path("results")
CATEGORY = "genome_mapped_but_transcriptome_unassigned"


def table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def save_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2), encoding="utf-8")


def marker(path):
    Path(path).write_text(time.strftime("%Y-%m-%dT%H:%M:%S%z") + "\n", encoding="ascii")


def event(stage, status, **details):
    (OUT / "logs").mkdir(parents=True, exist_ok=True)
    with (OUT / "logs/events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(dict(time=time.time(), stage=stage, status=status, **details)) + "\n")
    print(stage, status, details, flush=True)


def command(argv, name):
    argv = [str(x) for x in argv]
    env = command_environment(argv[0])
    event(name, "START", argv=argv, PATH_prefix=str(Path(argv[0]).parent))
    with (OUT / "logs" / (name + ".log")).open("w", encoding="utf-8") as log:
        subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
    event(name, "COMPLETE")


def command_environment(binary):
    env = os.environ.copy()
    env["PATH"] = str(Path(binary).parent) + os.pathsep + env.get("PATH", "")
    return env


def fragment_id(header):
    name = header.split()[0].lstrip(b"@")
    if name.endswith((b"/1", b"/2")):
        name = name[:-2]
    return name


def fastq_records(path):
    with gzip.open(path, "rb") as handle:
        while True:
            h = handle.readline()
            if not h:
                return
            s, p, q = handle.readline(), handle.readline(), handle.readline()
            if not h.startswith(b"@") or not p.startswith(b"+") or not q:
                raise ValueError("Malformed/truncated FASTQ: " + str(path))
            if len(s.rstrip()) != len(q.rstrip()) or not s.rstrip():
                raise ValueError("FASTQ sequence/quality mismatch: " + str(path))
            yield (h, s, p, q)


def paired_records(path1, path2):
    first, second = iter(fastq_records(path1)), iter(fastq_records(path2))
    while True:
        a, b = next(first, None), next(second, None)
        if a is None and b is None:
            return
        if a is None or b is None or fragment_id(a[0]) != fragment_id(b[0]):
            raise ValueError("FASTQ mate count/order/identity mismatch")
        yield fragment_id(a[0]), a, b


def expected_counts(config):
    rows = table(Path(config["revision"]) / "04_unmapped_reads/unmapped_read_categories.tsv")
    samples = {row["run"] for row in config["samples"]}
    values = {row["run"]: int(row["fragment_count"]) for row in rows if row["category"] == CATEGORY and row["run"] in samples}
    totals = [int(row["fragment_count"]) for row in rows if row["category"] == CATEGORY and row["run"] == "ALL_SAMPLES"]
    if len(values) != 9 or len(totals) != 1 or sum(values.values()) != totals[0]:
        raise ValueError("Category/sample counts do not reconcile")
    return values


def prepare(config, args):
    rev = Path(config["revision"])
    work = rev / "04_unmapped_reads"
    if os.path.commonpath([str(OUT.resolve()), str(rev.resolve())]) == str(rev.resolve()):
        raise ValueError("Outputs must not be inside the frozen revision")
    for name, binary in config["tools"].items():
        if not Path(binary).is_file() or not os.access(binary, os.X_OK):
            raise ValueError("Unavailable executable: " + name)
    version_args = {"getorf": ["-version"], "hmmsearch": ["-h"], "diamond": ["version"]}
    software = []
    for name, binary in config["tools"].items():
        process = subprocess.run([binary] + version_args.get(name, ["--version"]), capture_output=True, text=True, env=command_environment(binary))
        if process.returncode:
            raise ValueError("Executable version probe failed: " + name)
        software.append(dict(name=name, path=binary, version_output=(process.stdout+process.stderr)[:4000],
                             sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest()))
    save_json(OUT / "software_identity.json", software)
    expected = expected_counts(config)
    sources = [work / "unmapped_read_categories.tsv", work / "rnaseq_sample_sheet.tsv",
               work / "reference_annotation.gtf", Path(config["genome"]), Path(config["hmm_models"]),
               work / "classification_batch/fragment_categories.genome_mapped.tsv",
               work / "classification_batch/residual_1.fastq.gz", work / "classification_batch/residual_2.fastq.gz",
               work / "reference_databases/targeted_AMP_reference.faa",
               rev / "01_nslTP_curation/reference/negative_2S_seed_storage.faa",
               rev / "01_nslTP_curation/reference/ambiguous_prolamin_umbrella.faa",
               rev / "03_annotation_rescue/final_amp_primary_catalogue.tsv"]
    identities = []
    for path in sources:
        stat = path.stat()
        if stat.st_size == 0:
            raise ValueError("Empty input: " + str(path))
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if stat.st_size <= 40000000 else "NOT_REHASHED_LARGE_INPUT"
        identities.append(dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns, sha256=digest))
    for name in ["Genome", "SA", "SAindex", "genomeParameters.txt"]:
        if not (work / "star_index_sjdb149" / name).is_file():
            raise ValueError("STAR index incomplete")
    libraries = []
    for sample in config["samples"]:
        path = rev / "06_expression_reanalysis/salmon_quant" / sample["run"] / "lib_format_counts.json"
        data = json.loads(path.read_text())
        if data["expected_format"] != "ISR" or data["ISR"] <= 100 * max(1, data["ISF"]):
            raise ValueError("ISR assumption not supported: " + sample["run"])
        libraries.append(dict(run=sample["run"], ISR=data["ISR"], ISF=data["ISF"]))
    save_json(OUT / "source_identities.json", identities)
    save_json(OUT / "effective_config.json", config)
    save_json(OUT / "expected_fragment_counts.json", expected)
    save_json(OUT / "library_orientation.json", libraries)
    marker(OUT / "prepare.PASS")


def extract(config, args):
    batch = Path(config["revision"]) / "04_unmapped_reads/classification_batch"
    expected = expected_counts(config)
    states = {}
    category_counts = Counter()
    with (batch / "fragment_categories.genome_mapped.tsv").open("rb") as handle:
        for line in handle:
            name, category = line.rstrip(b"\r\n").split(b"\t")
            run = name.split(b".", 1)[0].decode("ascii")
            if category.decode() != CATEGORY or run not in expected or name in states:
                raise ValueError("Invalid/duplicate target category assignment")
            states[name] = False
            category_counts[run] += 1
    if dict(category_counts) != expected:
        raise ValueError("Target ID counts differ from frozen summary")
    (OUT / "reads").mkdir(exist_ok=True)
    handles = {(run, mate): gzip.open(OUT / "reads" / (run + "_" + str(mate) + ".fastq.gz"), "wb", compresslevel=1)
               for run in expected for mate in [1, 2]}
    selected, scanned = Counter(), 0
    try:
        for name, first, second in paired_records(batch / "residual_1.fastq.gz", batch / "residual_2.fastq.gz"):
            scanned += 1
            if name in states:
                if states[name]:
                    raise ValueError("Duplicate selected FASTQ fragment: " + name.decode())
                states[name] = True
                run = name.split(b".", 1)[0].decode("ascii")
                handles[run, 1].write(b"".join(first))
                handles[run, 2].write(b"".join(second))
                selected[run] += 1
            if scanned % 2000000 == 0:
                event("extract", "PROGRESS", scanned_pairs=scanned, selected_pairs=sum(selected.values()))
    finally:
        for handle in handles.values():
            handle.close()
    if dict(selected) != expected or not all(states.values()):
        raise ValueError("Incomplete selected fragment recovery")
    write_table(OUT / "fragment_extraction_audit.tsv", [dict(run=r, expected_fragments=expected[r], recovered_pairs=selected[r], status="PASS") for r in expected],
                ["run", "expected_fragments", "recovered_pairs", "status"])
    marker(OUT / "extract.PASS")


def align(config, args):
    tool = config["tools"]
    work = Path(config["revision"]) / "04_unmapped_reads"
    directory = OUT / "align" / args.run
    directory.mkdir(parents=True, exist_ok=True)
    command([tool["STAR"], "--runThreadN", args.threads, "--genomeDir", work / "star_index_sjdb149",
             "--genomeLoad", "NoSharedMemory", "--readFilesIn", OUT / "reads" / (args.run + "_1.fastq.gz"),
             OUT / "reads" / (args.run + "_2.fastq.gz"), "--readFilesCommand", "zcat", "--twopassMode", "Basic",
             "--outFileNamePrefix", str(directory) + "/", "--outSAMtype", "BAM", "SortedByCoordinate",
             "--limitBAMsortRAM", "16000000000", "--outSAMunmapped", "Within", "KeepPairs",
             "--outSAMattributes", "NH", "HI", "AS", "nM", "XS", "--outSAMstrandField", "intronMotif",
             "--outSAMattrRGline", "ID:" + args.run, "SM:" + args.run], args.run + ".STAR")
    bam = directory / "Aligned.sortedByCoord.out.bam"
    command([tool["samtools"], "quickcheck", "-v", bam], args.run + ".BAM_check")
    command([tool["samtools"], "index", "-@", max(1, args.threads - 1), bam], args.run + ".BAM_index")
    command([tool["samtools"], "flagstat", "-O", "json", bam], args.run + ".flagstat")
    marker(directory / "align.PASS")


def assemble(config, args):
    (OUT / "assembly").mkdir(exist_ok=True)
    command([config["tools"]["stringtie"], OUT / "align" / args.run / "Aligned.sortedByCoord.out.bam",
             "-G", Path(config["revision"]) / "04_unmapped_reads/reference_annotation.gtf", "--rf",
             "-p", args.threads, "-m", "90", "-c", "1", "-s", "2", "-f", "0.01", "-j", "2",
             "-a", "10", "-M", "0.5", "-l", "M6_" + args.run, "-o", OUT / "assembly" / (args.run + ".gtf")],
            args.run + ".StringTie")


def screen(config, args):
    from Bio import SeqIO
    tool, rev = config["tools"], Path(config["revision"])
    work = OUT / "screen"
    work.mkdir(exist_ok=True)
    paths = [str((OUT / "assembly" / (row["run"] + ".gtf")).resolve()) for row in config["samples"]]
    (work / "assemblies.txt").write_text("\n".join(paths) + "\n", encoding="utf-8")
    command([tool["stringtie"], "--merge", "-p", args.threads, "-m", "90", "-f", "0.01", "-l", "M6MERGE",
             "-o", work / "merged.gtf", work / "assemblies.txt"], "merge")
    command([tool["gffcompare"], "-r", rev / "04_unmapped_reads/reference_annotation.gtf", "-o", work / "comparison",
             work / "merged.gtf"], "gffcompare")
    command([tool["gffread"], work / "merged.gtf", "-g", config["genome"], "-w", work / "transcripts.fna"], "gffread")
    command([tool["getorf"], "-sequence", work / "transcripts.fna", "-outseq", work / "stop_delimited_orfs.faa",
             "-find", "0", "-minsize", "90", "-reverse", "N", "-table", "0", "-auto"], "getorf")
    proteins = list(SeqIO.parse(str(work / "stop_delimited_orfs.faa"), "fasta"))
    if not proteins:
        raise ValueError("No ORFs emitted: inspect assembly/ORF tool outputs before making a negative AMP claim")
    reference_files = [(rev / "04_unmapped_reads/reference_databases/targeted_AMP_reference.faa", ""),
                       (rev / "01_nslTP_curation/reference/negative_2S_seed_storage.faa", "CONTROL_2S|"),
                       (rev / "01_nslTP_curation/reference/ambiguous_prolamin_umbrella.faa", "CONTROL_PROLAMIN|")]
    seen, references = set(), []
    with (work / "labeled_references.faa").open("w", encoding="ascii") as handle:
        for source, prefix in reference_files:
            for rec in SeqIO.parse(str(source), "fasta"):
                rid = prefix + rec.id
                if rid in seen:
                    raise ValueError("Duplicate reference ID: " + rid)
                seen.add(rid)
                handle.write(">" + rid + "\n" + str(rec.seq) + "\n")
                references.append(dict(reference_id=rid, source=str(source), role="BOUNDARY_CONTROL" if prefix else "AMP_REFERENCE"))
    write_table(work / "reference_manifest.tsv", references, ["reference_id", "source", "role"])
    command([tool["diamond"], "makedb", "--in", work / "labeled_references.faa", "--db", work / "references"], "diamond_makedb")
    command([tool["diamond"], "blastp", "--query", work / "stop_delimited_orfs.faa", "--db", work / "references",
             "--out", work / "all_orf_reference_hits.tsv", "--sensitive", "--evalue", "0.001", "--id", "25",
             "--max-target-seqs", "0", "--threads", args.threads, "--outfmt", "6", "qseqid", "sseqid", "pident",
             "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "qcovhsp", "scovhsp", "evalue", "bitscore"], "diamond_ORFs")
    command([tool["hmmsearch"], "--cpu", max(1, args.threads - 1), "--cut_ga", "--noali", "--domtblout",
             work / "orf_family_domains.domtbl", config["hmm_models"], work / "stop_delimited_orfs.faa"], "hmmsearch_ORFs")
    marker(OUT / "screen.PASS")


def orf_coordinates(description):
    match = re.search(r"\[(\d+)\s*-\s*(\d+)\]", description)
    if match is None or int(match[1]) > int(match[2]):
        raise ValueError("Unexpected forward getorf header: " + description)
    return int(match[1]), int(match[2])


def overlaps(blocks, exons):
    return any(a < d and c < b for a, b in blocks for c, d in exons)


def support(bam, chrom, exons, strand):
    fragments, unique, multiple = set(), set(), set()
    covered = set()
    if chrom not in bam.references:
        raise ValueError("Chromosome missing from BAM: " + chrom)
    for read in bam.fetch(chrom, min(a for a, b in exons), max(b for a, b in exons)):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if not read.has_tag("NH"):
            raise ValueError("Missing NH tag in mapped alignment")
        rna_strand = "+" if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse) else "-"
        if strand in {"+", "-"} and rna_strand != strand:
            continue
        blocks = read.get_blocks()
        if not overlaps(blocks, exons):
            continue
        fragments.add(read.query_name)
        if read.get_tag("NH") == 1:
            unique.add(read.query_name)
            for a, b in blocks:
                for c, d in exons:
                    covered.update(range(max(a, c), min(b, d)))
        else:
            multiple.add(read.query_name)
    length = sum(b - a for a, b in exons)
    return dict(overlapping_fragments=len(fragments), unique_fragments=len(unique),
                multimapping_fragments=len(multiple), exon_bases=length, unique_covered_bases=len(covered),
                unique_exon_coverage_fraction=len(covered) / length if length else 0)


def evidence(config, args):
    import gffutils
    import pysam
    from Bio import SeqIO
    work = OUT / "screen"
    transcripts = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(str(work / "transcripts.fna"), "fasta")}
    proteins = {rec.id: rec for rec in SeqIO.parse(str(work / "stop_delimited_orfs.faa"), "fasta")}
    models, exons = {}, defaultdict(list)
    for f in gffutils.iterators.DataIterator(str(work / "comparison.annotated.gtf")):
        tid = f.attributes.get("transcript_id", [None])[0]
        if f.featuretype == "transcript":
            models[tid] = f
        elif f.featuretype == "exon":
            exons[tid].append((f.start - 1, f.end))
    domains = defaultdict(list)
    with (work / "orf_family_domains.domtbl").open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=22)
            domains[fields[0]].append(dict(model=fields[3], accession=fields[4], independent_evalue=fields[12], score=fields[13]))
    hits = defaultdict(list)
    with (work / "all_orf_reference_hits.tsv").open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            hits[fields[0]].append(fields)
    shortlist = sorted(set(domains) | {oid for oid, rows in hits.items() if any(not r[1].startswith("CONTROL_") for r in rows)})
    candidate_models = sorted({oid.rsplit("_", 1)[0] for oid in shortlist})
    if not set(candidate_models) <= set(models):
        raise ValueError("ORF/transcript/annotation identities do not join")
    evidence_rows, by_model = [], defaultdict(list)
    known_rows = []
    catalogue = table(Path(config["revision"]) / "03_annotation_rescue/final_amp_primary_catalogue.tsv")
    primary_gene_ids = {row["gene_id"] for row in catalogue}
    primary_transcript_ids = {row["protein_id"] for row in catalogue}
    tmap_paths = list(work.glob("comparison.*.tmap"))
    if len(tmap_paths) != 1:
        raise ValueError("Expected one gffcompare transcript relationship table")
    relationships = {row["qry_id"]: row for row in table(tmap_paths[0])}
    for sample in config["samples"]:
        run = sample["run"]
        with pysam.AlignmentFile(str(OUT / "align" / run / "Aligned.sortedByCoord.out.bam"), "rb") as bam:
            for tid in candidate_models:
                model = models[tid]
                stats = support(bam, model.seqid, sorted(exons[tid]), model.strand)
                row = dict(run=run, treatment=sample["treatment"], replicate=sample["replicate"], transcript_id=tid, **stats)
                evidence_rows.append(row)
                by_model[tid].append(row)
            for gene in catalogue:
                a, b = int(gene["start"]) - 1, int(gene["end"])
                chrom, strand = gene["chrom"], gene["strand"]
                for region, intervals in [("gene_body", [(a, b)]), ("flanks_10kb", [(max(0, a-10000), a), (b, min(b+10000, bam.get_reference_length(chrom)))])]:
                    intervals = [(x, y) for x, y in intervals if y > x]
                    stats = support(bam, chrom, intervals, strand)
                    renamed = {key.replace("exon_bases", "region_bases").replace("exon_coverage", "region_coverage"): value for key, value in stats.items()}
                    known_rows.append(dict(run=run, member_id=gene["member_id"], region=region, **renamed))
        event("evidence", "SAMPLE_COMPLETE", run=run)
    stats_fields = ["overlapping_fragments", "unique_fragments", "multimapping_fragments", "exon_bases", "unique_covered_bases", "unique_exon_coverage_fraction"]
    write_table(OUT / "candidate_per_sample_support.tsv", evidence_rows, ["run", "treatment", "replicate", "transcript_id"] + stats_fields)
    # For known gene-body/flank intervals the coverage columns describe those
    # intervals, not coding-exon coverage or proof of AMP transcription.
    region_fields = [key.replace("exon_bases", "region_bases").replace("exon_coverage", "region_coverage") for key in stats_fields]
    write_table(OUT / "known_AMP_region_support.tsv", known_rows, ["run", "member_id", "region"] + region_fields)
    result = []
    for oid in shortlist:
        tid = oid.rsplit("_", 1)[0]
        rec, seq, model = proteins[oid], transcripts[tid], models[tid]
        start, end = orf_coordinates(rec.description)
        protein = str(rec.seq)
        if end > len(seq) or end - start + 1 != len(protein) * 3:
            raise ValueError("ORF coordinates/length mismatch: " + oid)
        from Bio.Seq import Seq
        if str(Seq(seq[start-1:end]).translate()) != protein:
            raise ValueError("ORF translation does not match transcript: " + oid)
        positive = sorted([r for r in hits[oid] if not r[1].startswith("CONTROL_")], key=lambda r: -float(r[13]))
        negative = sorted([r for r in hits[oid] if r[1].startswith("CONTROL_")], key=lambda r: -float(r[13]))
        stops = seq[end:end+3] in {"TAA", "TAG", "TGA"}
        first_m = protein.find("M")
        complete_possible = stops and first_m >= 0 and len(protein) - first_m >= 30
        reliable_samples = sum(r["unique_fragments"] >= 3 and r["unique_exon_coverage_fraction"] >= 0.7 for r in by_model[tid])
        control_dominates = bool(negative) and (not positive or float(negative[0][13]) >= float(positive[0][13]))
        relation = relationships[tid]
        known_amp = relation["ref_gene_id"] in primary_gene_ids or relation["ref_id"] in primary_transcript_ids
        class_code = relation["class_code"]
        origin = "KNOWN_AMP_LOCUS" if known_amp else ("INTERGENIC_CANDIDATE" if class_code == "u" else "OTHER_OR_ALTERED_REFERENCE_RELATIONSHIP")
        priority = "KNOWN_AMP_REGION_SUPPORT" if known_amp else ("BOUNDARY_CONTROL_REVIEW" if control_dominates else ("FULL_READ_MODEL_VALIDATION_PRIORITY" if complete_possible and reliable_samples >= 2 else "FRAGMENT_OR_LOW_SUPPORT_REVIEW"))
        result.append(dict(orf_id=oid, transcript_id=tid, chrom=model.seqid, start=model.start, end=model.end, strand=model.strand,
            model_origin=origin, gffcompare_class=class_code, reference_gene_id=relation["ref_gene_id"],
            reference_id=relation["ref_id"], orf_start_in_transcript=start, orf_end_in_transcript=end,
            orf_length_aa=len(protein), downstream_stop_present=stops, first_methionine_position_aa=first_m+1 if first_m >= 0 else "NA",
            possible_start_stop_suborf=complete_possible, samples_with_3_unique_fragments_and_70pct_exon_coverage=reliable_samples,
            family_HMMs=";".join(d["accession"] for d in domains[oid]), best_AMP_reference=positive[0][1] if positive else "",
            best_AMP_bitscore=positive[0][13] if positive else "NA", best_boundary_control=negative[0][1] if negative else "",
            best_control_bitscore=negative[0][13] if negative else "NA", review_priority=priority,
            full_read_local_validation="NOT_YET_PERFORMED", final_gene_status="PROVISIONAL_NOT_ADDED_TO_CATALOGUE"))
    fields = list(result[0]) if result else ["orf_id", "transcript_id", "final_gene_status"]
    write_table(OUT / "candidate_evidence.tsv", result, fields)
    summary = dict(target_fragments=sum(expected_counts(config).values()), samples=len(config["samples"]),
        assembled_transcripts=len(transcripts), stop_delimited_orfs=len(proteins), HMM_or_AMP_similarity_candidates=len(result),
        candidate_transcript_models=len(candidate_models), review_priorities=dict(Counter(r["review_priority"] for r in result)),
        full_read_local_validation="NOT_YET_PERFORMED", catalogue_change="NONE", reviewer_M6_closed=False,
        scope="Genome-mapped, nondecoy, transcriptome-unassigned assembled transcript ORF screening; not exhaustive read-level or novel gene validation")
    save_json(OUT / "discovery_summary.json", summary)
    marker(OUT / "FOCUSED_DISCOVERY_COMPLETE.PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["prepare", "extract", "align", "assemble", "screen", "evidence"])
    parser.add_argument("--config", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--run")
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    OUT.mkdir(exist_ok=True)
    if hasattr(os, "nice"):
        os.nice(10)
    if args.run and args.run not in {r["run"] for r in config["samples"]}:
        raise ValueError("Unknown sample")
    event(args.stage, "START", run=args.run)
    try:
        globals()[args.stage](config, args)
    except Exception as exc:
        event(args.stage, "FAILED", run=args.run, error=str(exc))
        raise
    event(args.stage, "COMPLETE", run=args.run)


if __name__ == "__main__":
    main()
