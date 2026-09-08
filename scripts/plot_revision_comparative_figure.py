#!/usr/bin/env python3
"""Plot close-species AMP family counts and AMP-involving synteny."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath


FAMILIES = ["Defensin", "Snakin_GASA", "nsLTP"]
FAMILY_LABEL = {"Defensin": "Defensin", "Snakin_GASA": "Snakin/GASA", "nsLTP": "nsLTP"}
FAMILY_COLOR = {"Defensin": "#C43D3D", "Snakin_GASA": "#2E6F9E", "nsLTP": "#2F7D4A"}
SPECIES_COLOR = {"Reaumuria soongarica": "#4B4B4B", "Tamarix austromongolica": "#7A91A8"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--revision", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fai(path: Path) -> dict[str, int]:
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            name, length, *_ = line.rstrip("\n").split("\t")
            if name.startswith("Chr"):
                values[name] = int(length)
    return values


def chromosome_sort(name: str) -> tuple[int, str]:
    digits = "".join(character for character in name if character.isdigit())
    return (int(digits) if digits else 10**9, name)


def layout(lengths: dict[str, int], gap: float = 0.012) -> dict[str, tuple[float, float]]:
    chromosomes = sorted(lengths, key=chromosome_sort)
    usable = 1 - gap * (len(chromosomes) - 1)
    total = sum(lengths.values())
    positions = {}
    current = 0.0
    for chromosome in chromosomes:
        width = usable * lengths[chromosome] / total
        positions[chromosome] = (current, current + width)
        current += width + gap
    return positions


def genomic_x(chromosome: str, coordinate: float, lengths: dict[str, int], positions: dict[str, tuple[float, float]]) -> float:
    start, end = positions[chromosome]
    return start + (end - start) * coordinate / lengths[chromosome]


def curve(axis, x0: float, x1: float, y0: float, y1: float, color: str, alpha: float, linewidth: float) -> None:
    middle = (y0 + y1) / 2
    vertices = [(x0, y0), (x0, middle), (x1, middle), (x1, y1)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    axis.add_patch(patches.PathPatch(MplPath(vertices, codes), facecolor="none", edgecolor=color, alpha=alpha, linewidth=linewidth))


def main() -> None:
    args = parse_args()
    comp = args.revision / "05_comparative_genomics"
    counts = read_tsv(comp / "family_count_comparison.tsv")
    synteny = read_tsv(comp / "synteny_pairs.tsv")
    manifest = {row["combined_id"]: row for row in read_tsv(comp / "mcscanx_rso_tamarix/combined_id_manifest.tsv")}
    rso_lengths = fai(args.project / "data/genome.fa.fai")
    tau_lengths = fai(comp / "inputs/Tamarix_austromongolica/tau_genome.fasta.fai")
    rso_positions = layout(rso_lengths)
    tau_positions = layout(tau_lengths)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
        }
    )
    fig = plt.figure(figsize=(7.3, 9.2))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.34, 0.66], hspace=0.28)
    count_axis = fig.add_subplot(grid[0, 0])
    synteny_axis = fig.add_subplot(grid[1, 0])

    species = ["Reaumuria soongarica", "Tamarix austromongolica"]
    x = np.arange(len(FAMILIES))
    width = 0.34
    for offset, name in zip((-width / 2, width / 2), species):
        values = [int(next(row["high_confidence_count"] for row in counts if row["species"] == name and row["family"] == family)) for family in FAMILIES]
        bars = count_axis.bar(x + offset, values, width, color=SPECIES_COLOR[name], label=name)
        count_axis.bar_label(bars, fontsize=8, padding=2)
    count_axis.set_xticks(x)
    count_axis.set_xticklabels([FAMILY_LABEL[family] for family in FAMILIES], rotation=25, ha="right")
    count_axis.set_ylabel("Audited high-confidence members")
    count_axis.set_title("(A) Same-rule family counts", loc="left", fontweight="bold")
    count_axis.legend(frameon=False, fontsize=7.0)
    count_axis.spines[["top", "right"]].set_visible(False)
    count_axis.text(
        0,
        -0.24,
        "Tamarix defensin count includes five near-full-length homology models absent from its official annotation.",
        transform=count_axis.transAxes,
        fontsize=7.0,
        va="top",
        wrap=True,
    )

    bar_height = 0.055
    y_rso, y_tau = 0.08, 0.92
    for chromosome, (start, end) in rso_positions.items():
        synteny_axis.add_patch(patches.FancyBboxPatch((start, y_rso), end - start, bar_height, boxstyle="round,pad=0,rounding_size=0.006", facecolor="#4B4B4B", edgecolor="none", zorder=5))
        synteny_axis.text((start + end) / 2, y_rso - 0.018, chromosome, ha="center", va="top", fontsize=7.0)
    for chromosome, (start, end) in tau_positions.items():
        synteny_axis.add_patch(patches.FancyBboxPatch((start, y_tau - bar_height), end - start, bar_height, boxstyle="round,pad=0,rounding_size=0.006", facecolor="#7A91A8", edgecolor="none", zorder=5))
        synteny_axis.text((start + end) / 2, y_tau + 0.018, chromosome, ha="center", va="bottom", fontsize=7.0)

    direct = 0
    for row in synteny:
        rso = manifest.get(row["rso_gene"])
        tau = manifest.get(row["tamarix_gene"])
        if not rso or not tau or rso["chrom"] not in rso_positions or tau["chrom"] not in tau_positions:
            continue
        rso_x = genomic_x(rso["chrom"], (int(rso["start"]) + int(rso["end"])) / 2, rso_lengths, rso_positions)
        tau_x = genomic_x(tau["chrom"], (int(tau["start"]) + int(tau["end"])) / 2, tau_lengths, tau_positions)
        family = row["rso_amp_family"] or row["tamarix_amp_family"]
        is_direct = row["direct_amp_to_amp_anchor"] == "YES"
        direct += int(is_direct)
        curve(
            synteny_axis,
            rso_x,
            tau_x,
            y_rso + bar_height,
            y_tau - bar_height,
            FAMILY_COLOR.get(family, "#AFAFAF"),
            0.85 if is_direct else 0.22,
            1.25 if is_direct else 0.45,
        )
        if row["rso_amp_id"]:
            synteny_axis.plot(rso_x, y_rso + bar_height, marker="^", markersize=3.2, color=FAMILY_COLOR.get(family, "#777777"), zorder=6)
        if row["tamarix_amp_id"]:
            synteny_axis.plot(tau_x, y_tau - bar_height, marker="v", markersize=3.2, color=FAMILY_COLOR.get(family, "#777777"), zorder=6)

    synteny_axis.text(-0.02, y_rso + bar_height / 2, "R. soongarica", ha="right", va="center", fontstyle="italic", fontweight="bold")
    synteny_axis.text(-0.02, y_tau - bar_height / 2, "T. austromongolica", ha="right", va="center", fontstyle="italic", fontweight="bold")
    synteny_axis.set_xlim(-0.12, 1.01)
    synteny_axis.set_ylim(-0.03, 1.03)
    synteny_axis.axis("off")
    synteny_axis.set_title(f"(B) AMP-involving interspecies MCScanX anchors (direct AMP-AMP = {direct})", loc="left", fontweight="bold")
    handles = [patches.Patch(color=FAMILY_COLOR[family], label=FAMILY_LABEL[family]) for family in FAMILIES]
    handles.extend(
        [
            plt.Line2D([0], [0], color="#555555", linewidth=1.25, label="Direct AMP-AMP anchor"),
            plt.Line2D([0], [0], color="#AAAAAA", linewidth=0.45, alpha=0.5, label="AMP-to-non-AMP anchor in an interspecies block"),
        ]
    )
    synteny_axis.legend(handles=handles, frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.11), fontsize=7.0)

    fig.suptitle("Close-species AMP comparison in Tamaricaceae", fontsize=11, y=0.995)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    base = args.output_dir / "Figure8_Rso_Tamarix_comparative_genomics"
    for extension in ("svg", "pdf", "png"):
        fig.savefig(base.with_suffix(f".{extension}"))
    plt.close(fig)
    print(f"PASS comparative figure synteny_rows={len(synteny)} direct={direct}")


if __name__ == "__main__":
    main()
