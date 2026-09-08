# Promoter analysis update: 2026-09-08

This module supersedes the promoter occurrence/enrichment results in the original September 4/September 8 snapshot. Other analyses in the original package retain their own version and validation boundaries. The current input is the 63-gene primary catalogue, SHA-256 cddab7832f755995ca030594612d3097a77556b204c5aacb011eaac401f68e51; future catalogue changes require a rerun.

The analysis scans 2,000 bp upstream of annotated gene 5-prime boundaries, without claiming experimentally verified TSSs. The original 16 labels give 15 double-stranded patterns after merging TGACG/CGTCA. Exact sequence sources, custom-pattern status and functional-evidence limits are recorded in motif_sources.tsv and MOTIF_PROVENANCE.md. No new motifs were introduced.

| Scenario | Target genes | Background genes | Minimum BH-adjusted P | Enriched patterns |
| --- | ---: | ---: | ---: | ---: |
| Primary: complete 2-kb intervals | 63 | 21,698 | 0.790045 | 0 |
| S1: additionally A/C/G/T only and no other-gene overlap, in both cohorts | 60 | 19,959 | 0.874530 | 0 |
| S2: S1 plus annotated 5-prime UTR, in both cohorts | 51 | 19,146 | 0.756167 | 0 |

All 21,791 annotated-gene upstream sequences were reread from the genome. Primary target and background motif counts exactly reproduce the original 16-entry tables, before alias merging. An independent literal-string implementation verified 326,865 gene-pattern counts, all sequence hashes and all 45 statistical comparisons. Tests are one-sided Fisher greater; BH is reported within each 15-pattern scenario and additionally across all 45 comparisons. Complete results, including non-significant outcomes, are retained. Uniform QC is not full matching for GC content, expression or evolutionary dependence.

Download and extract [the complete promoter update](R1M8_promoter_revision_20260908.zip) for all input intervals, result tables, figures and checksums.

## Reproduce

Python >=3.8 standard library suffices for the analysis. Plotting additionally requires NumPy and Matplotlib (tested versions are recorded in figures/figure5_qa.json). Download the public reference genome/GFF using the original repository DATA_SOURCES.md and verify the input hashes in results/summary.json. The genome itself is not redistributed here.

```sh
python audit_promoters_v2.py --genome /path/to/genome.fa --gff /path/to/genome.gff --catalogue final_amp_primary_catalogue.tsv --outdir rerun_results
python verify_promoters_v2.py --results rerun_results --sources motif_sources.tsv --out independent_validation_rerun.json
python plot_promoters_v2.py --results rerun_results --outdir rerun_figures
```

Use a new output directory. Optional --legacy-results points to the original promoter result directory and asserts equality with its target/background counts. The included gzip FASTA contains all extracted upstream intervals in transcriptional orientation, allowing independent motif-count verification without redownloading the full genome. Individual coordinates, overlaps, ambiguity flags, UTR annotation flags and cohort membership are in results/all_gene_upstream_audit.tsv.

Figure 5 displays ln(1 + site count), with fixed family/member ordering, no clustering and no row scaling. The two annotation tracks mark other-gene overlap and absent 5-prime UTR annotation. A missing annotation is not biological absence of a UTR. Motif occurrence or enrichment does not establish binding, treatment-specific regulation or mechanism. No significant enrichment does not establish absence of regulation.

This update includes no manuscript, reviewer correspondence, credentials, raw RNA reads, or new wet-laboratory validation. Reuse permissions remain governed by LICENSE_NOTICE.md in the parent repository.
