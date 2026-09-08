import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = "/path/to/rsoamp/revision_R1_20260902";
const outDir = `${root}/09_manuscript_response/supplementary_workbook`;
const previewDir = `${outDir}/previews`;

const specs = [
  ["S1_Catalogue", "03_annotation_rescue/final_amp_primary_catalogue.tsv", "Primary 63-gene evidence-defined catalogue"],
  ["S2_AllAudited", "03_annotation_rescue/final_amp_catalogue_all_audited_members.tsv", "Primary, ambiguous, excluded, and genomic ORF candidates"],
  ["S3_nsLTP_Evidence", "01_nslTP_curation/nsLTP_evidence_matrix.tsv", "Candidate-level nsLTP/prolamin boundary evidence"],
  ["S4_nsLTP_Excluded", "01_nslTP_curation/excluded_candidates.tsv", "Audited nsLTP exclusions and reasons"],
  ["S5_BackfillLoci", "03_annotation_rescue/six_loci_audit.tsv", "Six historical backfill-locus audit"],
  ["S6_GenomeRescue", "03_annotation_rescue/genome_wide_rescue_hits.tsv", "Genome-wide rescue loci and annotation overlap"],
  ["S7_PhyloGroups", "02_phylogeny/rso_supported_group_sensitivity.tsv", "R. soongarica group recurrence across four strategies"],
  ["S8_PhyloSupport", "02_phylogeny/phylogeny_branch_support.tsv", "Branch-level SH-aLRT and UFBoot support"],
  ["S9_PromoterOccur", "07_promoter_kaks_methods/cis_element_full_results.tsv", "Promoter motif occurrence and coordinates"],
  ["S10_PromoterEnrich", "07_promoter_kaks_methods/cis_element_enrichment.tsv", "Fisher tests and BH-adjusted promoter enrichment"],
  ["S11_KaKsAudit", "07_promoter_kaks_methods/kaks_quality_audit.tsv", "All historical pairs with reliability gates"],
  ["S12_AMP_TPM", "06_expression_reanalysis/results/final_AMP_TPM.tsv", "Primary-catalogue transcript abundance (TPM)"],
  ["S13_AMP_NormCounts", "06_expression_reanalysis/results/final_AMP_normalized_counts.tsv", "DESeq2 normalized counts"],
  ["S14_AMP_DE", "06_expression_reanalysis/results/final_AMP_DE_results.tsv", "Complete primary-catalogue differential-expression results"],
  ["S15_DE_Summary", "06_expression_reanalysis/results/final_AMP_DE_summary.tsv", "Family-level DE counts by contrast"],
  ["S16_Mapping", "06_expression_reanalysis/salmon_quant_integrity.tsv", "Nine-library quantification integrity and mapping"],
  ["S17_UnmappedCat", "04_unmapped_reads/unmapped_read_categories.tsv", "Mutually exclusive unmapped-read categories"],
  ["S18_UnmappedAMP", "04_unmapped_reads/unmapped_amp_hits.tsv", "Targeted AMP-reference hits in residual reads"],
  ["S19_FamilyCounts", "05_comparative_genomics/family_count_comparison.tsv", "Evidence-matched family counts"],
  ["S20_Orthogroups", "05_comparative_genomics/orthogroups.tsv", "Reciprocal-best-hit groups involving AMP candidates"],
  ["S21_Synteny", "05_comparative_genomics/rso_tamarix_amp_synteny.tsv", "AMP-involving interspecies MCScanX anchors"],
  ["S22_qPCR_Candidates", "06_expression_reanalysis/qPCR_candidate_packet.tsv", "Wet-laboratory validation candidate packet"],
  ["S23_Mappability", "06_expression_reanalysis/results/key_gene_mappability.tsv", "Transcript mappability audit for prioritized genes"],
  ["S24_RefGenes", "06_expression_reanalysis/reference_gene_candidates.tsv", "Candidate qPCR reference genes; experimental stability testing required"],
  ["S25_Treatment", "06_expression_reanalysis/source_treatment_metadata.tsv", "Traceable source treatment metadata"],
];

function parseTsv(text) {
  const rows = text.replace(/^\uFEFF/, "").replace(/\r/g, "").split("\n").filter((x) => x.length > 0);
  return rows.map((line) => line.split("\t").map(coerce));
}

function coerce(value) {
  if (value === "") return null;
  if (/^(?:[-+]?\d+(?:\.\d+)?|[-+]?\.\d+)(?:[eE][-+]?\d+)?$/.test(value) && !/^0\d+$/.test(value)) {
    const n = Number(value);
    if (Number.isFinite(n)) return n;
  }
  return value;
}

function colLetter(n) {
  let s = "";
  while (n > 0) {
    n -= 1;
    s = String.fromCharCode(65 + (n % 26)) + s;
    n = Math.floor(n / 26);
  }
  return s;
}

function widthFor(values) {
  const maxLen = Math.max(...values.map((v) => String(v ?? "").length));
  return Math.max(10, Math.min(34, maxLen + 2));
}

function safeTableName(sheetName) {
  return `Tbl_${sheetName.replace(/[^A-Za-z0-9_]/g, "_")}`.slice(0, 240);
}

const workbook = Workbook.create();
const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
const readmeRows = [
  ["Genes RsoAMP major-revision supplementary tables", "Value / evidence boundary"],
  ["Generated", "2026-09-02"],
  ["Primary catalogue", "63 complete annotated genes: 9 defensins, 18 Snakin/GASA, 36 canonical nsLTPs"],
  ["Historical backfills", "6 genomic AMP-like ORF candidates; not counted as genes and excluded from expression, promoter, duplication, and synteny analyses"],
  ["nsLTP boundary", "RsLTP17 and RsLTP19 retained as predicted LTPg; RsLTP2, RsLTP18, RsLTP20, and RsLTP30 excluded from the canonical catalogue"],
  ["Phylogeny", "Sensitivity-qualified supported clades only; no complete formal subfamily classification"],
  ["Expression threshold", "BH-adjusted P < 0.05 and absolute unshrunken MLE log2 fold change > 1"],
  ["Treatment", "0, 200, and 400 mM Na2SO4; three biological replicates per condition"],
  ["qRT-PCR", "Not performed. Candidate and reference-gene sheets are validation planning outputs, not experimental confirmation"],
  ["Public repository", "File S1 is repository-ready; no permanent public URL has been assigned"],
  ["Genome source", "https://doi.org/10.1038/s41597-024-03644-y ; https://doi.org/10.6084/m9.figshare.25533064.v2"],
  ["RNA-seq source", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1063761"],
  ["Comparator source", "https://doi.org/10.1093/dnares/dsae021 ; https://doi.org/10.6084/m9.figshare.25106726"],
  ["Sheet", "Source file and purpose"],
  ...specs.map(([sheet, rel, purpose]) => [sheet, `${rel} | ${purpose}`]),
];
readme.getRangeByIndexes(0, 0, readmeRows.length, 2).values = readmeRows;
readme.getRange("A1:B1").format = { fill: "#1F4E78", font: { bold: true, color: "#FFFFFF", size: 14 }, wrapText: true };
readme.getRange("A14:B14").format = { fill: "#548235", font: { bold: true, color: "#FFFFFF" } };
readme.getRangeByIndexes(1, 0, readmeRows.length - 1, 2).format = { font: { name: "Arial", size: 10 }, verticalAlignment: "top", wrapText: true };
readme.getRangeByIndexes(0, 0, readmeRows.length, 2).format.borders = { preset: "inside", style: "thin", color: "#D9E2F3" };
readme.getRange("A:A").format.columnWidth = 25;
readme.getRange("B:B").format.columnWidth = 95;
readme.freezePanes.freezeRows(1);

for (const [sheetName, rel, purpose] of specs) {
  const sheet = workbook.worksheets.add(sheetName);
  sheet.showGridLines = false;
  const file = path.join(root, rel);
  let rows;
  try {
    rows = parseTsv(await fs.readFile(file, "utf8"));
  } catch (error) {
    rows = [["status", "source_file", "purpose", "note"], ["PENDING", rel, purpose, "Source output was not available when the workbook was generated; regenerate before submission."]];
  }
  if (rows.length === 1) rows.push(Array(rows[0].length).fill(null));
  const colCount = Math.max(...rows.map((r) => r.length));
  rows = rows.map((r) => r.concat(Array(colCount - r.length).fill(null)));
  sheet.getRangeByIndexes(0, 0, rows.length, colCount).values = rows;
  const used = sheet.getRangeByIndexes(0, 0, rows.length, colCount);
  used.format = { font: { name: "Arial", size: 9 }, verticalAlignment: "top" };
  const header = sheet.getRangeByIndexes(0, 0, 1, colCount);
  header.format = {
    fill: "#1F4E78",
    font: { bold: true, color: "#FFFFFF", size: 9 },
    wrapText: true,
    verticalAlignment: "center",
    borders: { preset: "outside", style: "medium", color: "#1F4E78" },
  };
  header.format.rowHeight = 32;
  if (rows.length > 1) {
    const body = sheet.getRangeByIndexes(1, 0, rows.length - 1, colCount);
    body.format.borders = { preset: "inside", style: "thin", color: "#E7E6E6" };
    const table = sheet.tables.add(`A1:${colLetter(colCount)}${rows.length}`, true, safeTableName(sheetName));
    table.style = "TableStyleMedium2";
    table.showBandedRows = true;
  }
  for (let c = 0; c < colCount; c += 1) {
    const values = rows.slice(0, Math.min(rows.length, 250)).map((r) => r[c]);
    sheet.getRangeByIndexes(0, c, rows.length, 1).format.columnWidth = widthFor(values);
  }
  sheet.freezePanes.freezeRows(1);
}

await fs.mkdir(previewDir, { recursive: true });
for (const [index, sheet] of workbook.worksheets.items.entries()) {
  const used = sheet.getUsedRange(true);
  const maxRows = Math.min(used.rowCount, 45);
  const maxCols = Math.min(used.columnCount, 16);
  const range = `A1:${colLetter(maxCols)}${maxRows}`;
  const png = await workbook.render({ sheetName: sheet.name, range, scale: 0.85, format: "png" });
  await fs.writeFile(`${previewDir}/${String(index + 1).padStart(2, "0")}_${sheet.name}.png`, new Uint8Array(await png.arrayBuffer()));
}

const summary = await workbook.inspect({ kind: "sheet", include: "id,name", maxChars: 12000 });
console.log(summary.ndjson);
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 100 }, summary: "final formula error scan" });
console.log(errors.ndjson);

await fs.mkdir(outDir, { recursive: true });
const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(`${outDir}/Genes_RsoAMP_Supplementary_Tables_R1_20260902.xlsx`);
