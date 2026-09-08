#!/usr/bin/env python3
"""Generate editable revision figures 1, 3, 4, 5, 6 and 7 from frozen tables."""

from __future__ import annotations

import argparse
import csv
import math
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


FAMILY_ORDER = ["Defensin", "Snakin_GASA", "nsLTP"]
FAMILY_LABEL = {"Defensin": "Defensin", "Snakin_GASA": "Snakin/GASA", "nsLTP": "nsLTP"}
FAMILY_COLOR = {"Defensin": "#C43D3D", "Snakin_GASA": "#2E6F9E", "nsLTP": "#2F7D4A"}
CONDITION_COLOR = {"CK": "#4C4C4C", "S200": "#D8903A", "S400": "#B23A48"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--revision", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def configure() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
        }
    )


def rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def finite_number(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def member_number(member_id: str) -> int:
    match = re.search(r"(\d+)$", member_id)
    return int(match.group(1)) if match else 10**9


def member_sort(row: dict[str, str]) -> tuple[int, int, str]:
    return (FAMILY_ORDER.index(row["family"]), member_number(row["member_id"]), row["member_id"])


def save(fig, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"))
    plt.close(fig)


def plot_chromosomes(project: Path, catalogue: list[dict[str, str]], output: Path) -> None:
    lengths = {}
    with (project / "data/genome.fa.fai").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if re.fullmatch(r"Chr\d+", fields[0]):
                lengths[fields[0]] = int(fields[1])
    chromosomes = sorted(lengths, key=lambda value: int(value[3:] if value.startswith("Chr") else value))
    maximum_mb = max(lengths.values()) / 1e6
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 10.4), sharey=True)
    chromosome_groups = [chromosomes[:6], chromosomes[6:]]
    half_width = 0.12
    minimum_gap = maximum_mb * 0.018
    for panel, (axis, group) in enumerate(zip(axes, chromosome_groups), start=1):
        for index, chromosome in enumerate(group):
            chromosome_length = lengths[chromosome] / 1e6
            axis.add_patch(
                patches.FancyBboxPatch(
                    (index - half_width, 0), 2 * half_width, chromosome_length,
                    boxstyle="round,pad=0.015,rounding_size=0.14",
                    linewidth=0.7, edgecolor="#555555", facecolor="#E7E7E7",
                )
            )
            axis.text(index, -3, chromosome, ha="center", va="top", fontsize=8, fontweight="bold")
            members = sorted(
                [row for row in catalogue if row["chrom"] == chromosome],
                key=lambda row: int(row["start"]),
            )
            genomic_positions = [int(row["start"]) / 1e6 for row in members]
            label_positions = list(genomic_positions)
            for position in range(1, len(label_positions)):
                label_positions[position] = max(label_positions[position], label_positions[position - 1] + minimum_gap)
            if label_positions and label_positions[-1] > chromosome_length:
                shift = label_positions[-1] - chromosome_length
                label_positions = [value - shift for value in label_positions]
                for position in range(len(label_positions) - 2, -1, -1):
                    label_positions[position] = min(
                        label_positions[position], label_positions[position + 1] - minimum_gap
                    )
            if label_positions and label_positions[0] < 0:
                shift = -label_positions[0]
                label_positions = [value + shift for value in label_positions]
            for member, genomic, label in zip(members, genomic_positions, label_positions):
                color = FAMILY_COLOR[member["family"]]
                axis.plot([index - half_width, index + half_width], [genomic, genomic], color=color, linewidth=1.2)
                axis.plot([index + half_width, index + 0.20], [genomic, label], color=color, linewidth=0.45)
                axis.text(
                    index + 0.22, label, member["member_id"],
                    ha="left", va="center", fontsize=7.0, color=color,
                )
        axis.set_ylim(maximum_mb + 5, -9)
        axis.set_xlim(-0.55, len(group) - 0.15)
        axis.set_xticks([])
        axis.spines[["top", "right", "bottom"]].set_visible(False)
        axis.set_title(
            f"({'A' if panel == 1 else 'B'}) {group[0]}-{group[-1]}",
            loc="left", fontsize=9, fontweight="bold",
        )
    fig.supylabel("Position (Mb)", x=0.01, fontsize=8)
    handles = [patches.Patch(color=FAMILY_COLOR[family], label=FAMILY_LABEL[family]) for family in FAMILY_ORDER]
    axes[-1].legend(handles=handles, loc="lower left", frameon=False, ncol=3, title="Primary catalogue")
    fig.suptitle(f"Chromosomal distribution of audited AMP genes (n = {len(catalogue)})", y=0.995)
    fig.tight_layout(rect=(0.02, 0.01, 1, 0.98))
    save(fig, output / "Figure1_chromosomal_distribution")


def parse_cds(gff: Path) -> dict[str, list[tuple[int, int]]]:
    cds = defaultdict(list)
    with gff.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].strip(";").split(";") if "=" in item)
            for parent in attrs.get("Parent", "").split(","):
                if parent:
                    cds[parent].append((int(fields[3]), int(fields[4])))
    return cds


def plot_structure(project: Path, catalogue: list[dict[str, str]], output: Path) -> None:
    cds = parse_cds(project / "data/genome.gff")
    fig, axis = plt.subplots(figsize=(7.3, max(8.0, 0.20 * len(catalogue))))
    for y, member in enumerate(catalogue):
        segments = sorted(cds.get(member["protein_id"], []))
        if not segments:
            continue
        minimum = min(start for start, _ in segments)
        maximum = max(end for _, end in segments)
        span = max(1, maximum - minimum)
        axis.plot([0, 1], [y, y], color="#9A9A9A", linewidth=0.55)
        for start, end in segments:
            axis.add_patch(
                patches.Rectangle(
                    ((start - minimum) / span, y - 0.3),
                    max((end - start + 1) / span, 0.006),
                    0.6,
                    facecolor=FAMILY_COLOR[member["family"]], edgecolor="none",
                )
            )
    axis.set_yticks(range(len(catalogue)))
    axis.set_yticklabels([row["member_id"] for row in catalogue], fontsize=8.0)
    axis.invert_yaxis()
    axis.set_xlim(-0.02, 1.02)
    axis.set_xlabel("Normalized gene span (boxes: coding exons; lines: introns)")
    axis.set_title("Gene structures of primary-catalogue AMP members")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(left=False)
    save(fig, output / "Figure3_gene_structure")


def plot_motifs(project: Path, catalogue: list[dict[str, str]], output: Path) -> None:
    panel_heights = [
        max(2.4, 0.22 * sum(row["family"] == family for row in catalogue))
        for family in FAMILY_ORDER
    ]
    fig, axes = plt.subplots(
        3, 1,
        figsize=(7.3, sum(panel_heights) + 1.2),
        gridspec_kw={"height_ratios": panel_heights, "hspace": 0.42},
    )
    palette = plt.cm.tab10.colors
    for axis, family in zip(axes, FAMILY_ORDER):
        xml_path = project / f"intermediate/meme/{family}/meme.xml"
        root = ET.parse(xml_path).getroot()
        sequences = {entry.get("id"): (entry.get("name"), int(entry.get("length"))) for entry in root.iter("sequence")}
        motif_defs = {entry.get("id"): (entry.get("name") or entry.get("id"), int(entry.get("width"))) for entry in root.find("motifs")}
        motif_order = list(motif_defs)
        sites = defaultdict(list)
        for scanned in root.iter("scanned_sites"):
            for site in scanned.findall("scanned_site"):
                sites[scanned.get("sequence_id")].append((site.get("motif_id"), int(site.get("position"))))
        id_by_name = {name: identifier for identifier, (name, _) in sequences.items()}
        family_members = [row for row in catalogue if row["family"] == family]
        for y, member in enumerate(family_members):
            identifier = id_by_name.get(member["member_id"])
            if identifier is None:
                continue
            length = sequences[identifier][1]
            axis.plot([0, length], [y, y], color="#C5C5C5", linewidth=0.6)
            for motif, position in sites.get(identifier, []):
                width = motif_defs[motif][1]
                axis.add_patch(
                    patches.Rectangle((position, y - 0.32), width, 0.64, facecolor=palette[motif_order.index(motif) % 10], edgecolor="none")
                )
        axis.set_yticks(range(len(family_members)))
        axis.set_yticklabels([row["member_id"] for row in family_members], fontsize=8.0)
        axis.invert_yaxis()
        axis.set_title(FAMILY_LABEL[family], color=FAMILY_COLOR[family])
        axis.set_xlabel("Precursor position (aa)")
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(left=False)
    motif_handles = [patches.Patch(color=palette[index], label=f"Motif {index + 1}") for index in range(5)]
    fig.suptitle("Conserved MEME motifs in primary-catalogue AMP proteins", y=0.998)
    fig.legend(handles=motif_handles, loc="upper center", bbox_to_anchor=(0.5, 0.976), ncol=5, frameon=False)
    fig.subplots_adjust(top=0.94)
    save(fig, output / "Figure4_conserved_motifs")


def plot_promoters(revision: Path, catalogue: list[dict[str, str]], output: Path) -> None:
    promoter_rows = rows(revision / "07_promoter_kaks_methods/cis_element_full_results.tsv")
    elements = []
    values = defaultdict(dict)
    for row in promoter_rows:
        if row["element"] not in elements:
            elements.append(row["element"])
        values[row["member_id"]][row["element"]] = float(row["nonredundant_genomic_site_count"])
    matrix = np.array([[values[row["member_id"]].get(element, 0) for element in elements] for row in catalogue], dtype=float)
    fig, axis = plt.subplots(figsize=(7.3, max(8.0, 0.20 * len(catalogue))))
    image = axis.pcolormesh(
        np.arange(len(elements) + 1),
        np.arange(len(catalogue) + 1),
        np.log1p(matrix),
        cmap="YlOrRd",
        shading="flat",
    )
    axis.set_xlim(0, len(elements))
    axis.set_ylim(len(catalogue), 0)
    axis.set_xticks(np.arange(len(elements)) + 0.5)
    axis.set_xticklabels(elements, rotation=50, ha="right", fontsize=8.0)
    axis.set_yticks(np.arange(len(catalogue)) + 0.5)
    axis.set_yticklabels([row["member_id"] for row in catalogue], fontsize=8.0)
    colorbar = fig.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
    if colorbar.solids is not None:
        colorbar.solids.set_rasterized(False)
    colorbar.set_label("log(1 + nonredundant site count)")
    axis.set_title("Sequence-level cis-element occurrences in 2-kb upstream regions")
    save(fig, output / "Figure5_promoter_motif_occurrence")


def plot_expression(revision: Path, catalogue: list[dict[str, str]], output: Path) -> None:
    matrix_rows = rows(revision / "06_expression_reanalysis/results/final_AMP_VST.tsv")
    sample_order = ["SRR27540881", "SRR27540880", "SRR27540879", "SRR27540878", "SRR27540877", "SRR27540876", "SRR27540875", "SRR27540884", "SRR27540883"]
    condition = {sample: ("CK" if index < 3 else "S200" if index < 6 else "S400") for index, sample in enumerate(sample_order)}
    labels, families, matrix = [], [], []
    for row in matrix_rows:
        if row["quantification_status"] != "QUANTIFIED":
            continue
        labels.append(row["member_id"])
        families.append(row["family"])
        matrix.append([float(row[sample]) for sample in sample_order])
    matrix = np.array(matrix, dtype=float)
    means = matrix.mean(axis=1, keepdims=True)
    standard = matrix.std(axis=1, keepdims=True)
    standard[standard == 0] = 1
    zscores = (matrix - means) / standard
    figure_height = max(8.0, 0.20 * len(labels))
    fig, (family_axis, heat_axis) = plt.subplots(1, 2, figsize=(7.3, figure_height), gridspec_kw={"width_ratios": [0.04, 1], "wspace": 0.02})
    for index, family in enumerate(families):
        family_axis.add_patch(patches.Rectangle((0, len(labels) - index - 1), 1, 1, color=FAMILY_COLOR[family], linewidth=0))
    family_axis.set_xlim(0, 1)
    family_axis.set_ylim(0, len(labels))
    family_axis.axis("off")
    limit = np.nanpercentile(np.abs(zscores), 98) or 2
    image = heat_axis.pcolormesh(
        np.arange(zscores.shape[1] + 1),
        np.arange(zscores.shape[0] + 1),
        zscores,
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
        shading="flat",
    )
    heat_axis.set_xlim(0, zscores.shape[1])
    heat_axis.set_ylim(zscores.shape[0], 0)
    heat_axis.set_yticks(np.arange(len(labels)) + 0.5)
    heat_axis.set_yticklabels(labels, fontsize=8.0)
    heat_axis.set_xticks(np.arange(9) + 0.5)
    heat_axis.set_xticklabels(["rep1", "rep2", "rep3"] * 3, rotation=90)
    heat_axis.axvline(3, color="white", linewidth=1.4)
    heat_axis.axvline(6, color="white", linewidth=1.4)
    for x, name in ((1.5, "CK"), (4.5, "S200"), (7.5, "S400")):
        heat_axis.text(
            x,
            1.006,
            name,
            transform=heat_axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            color=CONDITION_COLOR[name],
            fontweight="bold",
            clip_on=False,
        )
    heat_axis.tick_params(length=0)
    colorbar = fig.colorbar(image, ax=heat_axis, fraction=0.026, pad=0.02)
    if colorbar.solids is not None:
        colorbar.solids.set_rasterized(False)
    colorbar.set_label("Row z-score of VST counts")
    heat_axis.set_title(r"AMP expression across the Na$_2$SO$_4$ treatment series", pad=28)
    save(fig, output / "Figure6_expression_heatmap")

    de_rows = rows(revision / "06_expression_reanalysis/results/final_AMP_DE_results.tsv")
    contrasts = ["S200_vs_CK", "S400_vs_CK"]
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 4.6), sharey=True)
    for axis, contrast in zip(axes, contrasts):
        selected = [
            row
            for row in de_rows
            if row["contrast"] == contrast
            and row["quantification_status"] == "QUANTIFIED"
            and finite_number(row["log2FoldChange"])
            and finite_number(row["padj_BH"])
        ]
        significant = []
        for row in selected:
            fold = float(row["log2FoldChange"])
            adjusted = float(row["padj_BH"])
            is_de = row["DE_call_FDR_lt_0.05_abs_log2FC_gt_1"].upper() == "TRUE"
            color = FAMILY_COLOR[row["family"]] if is_de else "#CFCFCF"
            axis.scatter(fold, -math.log10(max(adjusted, 1e-300)), s=18 if is_de else 9, color=color, alpha=0.9 if is_de else 0.55, edgecolors="none")
            if is_de:
                significant.append(row)
        for index, row in enumerate(sorted(significant, key=lambda value: float(value["padj_BH"]))[:10]):
            x_offset = -4 if index % 2 == 0 else 4
            axis.annotate(
                row["member_id"],
                (float(row["log2FoldChange"]), -math.log10(max(float(row["padj_BH"]), 1e-300))),
                xytext=(x_offset, 5 + 3 * (index % 3)),
                textcoords="offset points",
                fontsize=7.0,
                ha="right" if x_offset < 0 else "left",
                va="bottom",
            )
        axis.axhline(-math.log10(0.05), color="#777777", linestyle="--", linewidth=0.6)
        axis.axvline(-1, color="#777777", linestyle="--", linewidth=0.6)
        axis.axvline(1, color="#777777", linestyle="--", linewidth=0.6)
        axis.set_xlabel("log2 fold change (MLE)")
        axis.set_title(f"{contrast.replace('_vs_', ' vs ')} (DE = {len(significant)})")
    axes[0].set_ylabel("-log10 Benjamini-Hochberg adjusted P")
    handles = [patches.Patch(color=FAMILY_COLOR[family], label=FAMILY_LABEL[family]) for family in FAMILY_ORDER]
    handles.append(patches.Patch(color="#CFCFCF", label="Not significant"))
    axes[1].legend(handles=handles, frameon=False, loc="upper right")
    fig.suptitle("Differential expression of primary-catalogue AMP genes")
    fig.tight_layout()
    save(fig, output / "Figure7_expression_volcano")


def main() -> None:
    args = parse_args()
    configure()
    catalogue = sorted(rows(args.revision / "03_annotation_rescue/final_amp_primary_catalogue.tsv"), key=member_sort)
    plot_chromosomes(args.project, catalogue, args.output_dir)
    plot_structure(args.project, catalogue, args.output_dir)
    plot_motifs(args.project, catalogue, args.output_dir)
    plot_promoters(args.revision, catalogue, args.output_dir)
    plot_expression(args.revision, catalogue, args.output_dir)
    print(f"PASS core revision figures primary_catalogue={len(catalogue)}")


if __name__ == "__main__":
    main()
