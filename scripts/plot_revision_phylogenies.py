#!/usr/bin/env python3
"""Plot large editable representative trees without changing inferred topology."""

from __future__ import annotations

import argparse
import copy
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from Bio import Phylo


FAMILIES = ["Defensin", "Snakin_GASA", "nsLTP"]
LABELS = {"Defensin": "Defensin", "Snakin_GASA": "Snakin/GASA", "nsLTP": "nsLTP"}
COLORS = {
    "primary_Rso": "#B33A3A",
    "excluded_Rso": "#D47A2C",
    "Tamarix": "#2E6F9E",
    "positive_reference": "#555555",
    "boundary_control": "#8C6D9E",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--strategy", default="full.gappyout")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def support_label(text: str | None) -> str | None:
    if not text:
        return None
    match = re.fullmatch(r"([0-9.]+)/([0-9.]+)", str(text))
    if not match:
        return None
    sh, uf = map(float, match.groups())
    return f"{sh:.0f}/{uf:.0f}" if sh >= 80 and uf >= 95 else None


def tip_class(name: str, primary_members: set[str], current_rso: set[str]) -> str:
    if name.startswith("RSO_"):
        label = name.removeprefix("RSO_")
        return "primary_Rso" if label in primary_members else "excluded_Rso"
    if name.startswith("TAU_"):
        return "Tamarix"
    if name.startswith(("NEG2S_", "AMBPRO_")):
        return "boundary_control"
    return "positive_reference"


def display_label(name: str, primary_members: set[str], current_rso: set[str]) -> str:
    if name.startswith("RSO_"):
        label = name.removeprefix("RSO_")
        if label in primary_members:
            return label
        if label in current_rso:
            return f"{label} [excluded]"
        return ""
    if name.startswith("TAU_"):
        return name.removeprefix("TAU_")
    if name.startswith("UNI_"):
        return name.removeprefix("UNI_")
    if name.startswith("NSPOS_"):
        return name.removeprefix("NSPOS_")
    if name.startswith("NEG2S_"):
        return name.removeprefix("NEG2S_") + " [2S control]"
    if name.startswith("AMBPRO_"):
        return ""
    return name


def main_display_tree(tree, retained_rso_labels: set[str]):
    """Retain primary/ambiguous RSO tips plus each tip's nearest context."""
    rso_tips = [
        tip
        for tip in tree.get_terminals()
        if (tip.name or "").startswith("RSO_")
        and (tip.name or "").removeprefix("RSO_") in retained_rso_labels
    ]
    context_tips = [
        tip
        for tip in tree.get_terminals()
        if not (tip.name or "").startswith("RSO_")
        and not (tip.name or "").startswith("AMBPRO_")
    ]
    keep = {tip.name for tip in rso_tips}
    for tip in rso_tips:
        nearest = min(context_tips, key=lambda candidate: (tree.distance(tip, candidate), candidate.name or ""))
        keep.add(nearest.name)
    display = copy.deepcopy(tree)
    for tip in list(display.get_terminals()):
        if tip.name not in keep:
            display.prune(tip)
    root = display.root
    for clade in list(display.get_nonterminals(order="postorder")):
        if clade is root:
            continue
        if support_label(clade.name) is None:
            display.collapse(clade)
    display.ladderize()
    return display, keep


def draw_tree(tree, base: Path, title: str, primary_members: set[str], current_rso: set[str], full: bool) -> None:
    rendered = {
        tip.name: display_label(tip.name or "", primary_members, current_rso)
        for tip in tree.get_terminals()
    }
    label_to_tip = {label: name for name, label in rendered.items() if label}

    def label_func(clade):
        if clade.is_terminal():
            return rendered.get(clade.name or "", "")
        return support_label(clade.name) if full else None

    def label_color(label):
        name = label_to_tip.get(label)
        if name is None:
            return "#333333"
        return COLORS[tip_class(name, primary_members, current_rso)]

    n_tips = tree.count_terminals()
    if full:
        width = 10.5
        height = max(7.0, 0.105 * n_tips)
        label_size = 5.7 if n_tips > 100 else 6.2
    else:
        width = 7.3
        height = max(7.0, 0.16 * n_tips)
        label_size = 8.0
    fig, axis = plt.subplots(figsize=(width, height))
    Phylo.draw(
        tree,
        axes=axis,
        do_show=False,
        label_func=label_func,
        label_colors=label_color,
        show_confidence=False,
    )
    for text in axis.texts:
        if re.fullmatch(r"\s*\d+/\d+\s*", text.get_text()):
            text.set_fontsize(4.8 if not full else 4.5)
            text.set_bbox({"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.15})
        else:
            text.set_fontsize(label_size)
    axis.set_title(title, loc="left", fontweight="bold")
    axis.set_xlabel("Substitutions per site")
    axis.set_ylabel("")
    axis.set_yticks([])
    axis.spines[["top", "right", "left"]].set_visible(False)
    legend_labels = {
        "primary_Rso": "R. soongarica primary catalogue",
        "excluded_Rso": "R. soongarica audited exclusion",
        "Tamarix": "Tamarix austromongolica",
        "positive_reference": "Reviewed family reference",
        "boundary_control": "Explicit prolamin/2S boundary control",
    }
    present_classes = {
        tip_class(name, primary_members, current_rso)
        for name, label in rendered.items()
        if label
    }
    legend = [
        patches.Patch(color=COLORS[category], label=legend_labels[category])
        for category in legend_labels
        if category in present_classes
    ]
    fig.subplots_adjust(bottom=0.20 if not full else 0.15)
    axis.legend(
        handles=legend,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(0, -0.075),
        ncol=2 if not full else 3,
        fontsize=7.0 if not full else 6.0,
    )
    footer = (
        "Full inferred topology; node labels shown only for SH-aLRT >= 80 and UFBoot >= 95."
        if full
        else "Primary and named ambiguous R. soongarica candidates are shown with nearest context; all retained bifurcations meet SH-aLRT >= 80 and UFBoot >= 95, and unsupported nodes are collapsed."
    )
    fig.text(0.5, 0.002, footer, ha="center", va="bottom", fontsize=7.0 if not full else 6.5)
    for extension in ("svg", "pdf", "png"):
        fig.savefig(base.with_suffix(f".{extension}"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    primary = read_tsv(args.revision / "03_annotation_rescue/final_amp_primary_catalogue.tsv")
    primary_by_family = {
        family: {row["member_id"] for row in primary if row["family"] == family}
        for family in FAMILIES
    }
    reclassified = read_tsv(args.revision / "01_nslTP_curation/reclassified_candidates.tsv")
    named_ambiguous_nsltp = {
        row["member_id"]
        for row in reclassified
        if row.get("member_id")
        and row.get("final_classification") == "B_AMBIGUOUS_NSLTP_PROLAMIN_BOUNDARY"
    }
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
        }
    )
    manifest_rows = []
    for panel, family in zip("ABC", FAMILIES):
        tree_path = args.revision / f"02_phylogeny/{family}/{family}.{args.strategy}.treefile"
        full_tree = Phylo.read(tree_path, "newick")
        full_tree.ladderize()
        primary_members = primary_by_family[family]
        current_rso = {
            (tip.name or "").removeprefix("RSO_")
            for tip in full_tree.get_terminals()
            if (tip.name or "").startswith("RSO_")
        }
        retained_rso = set(primary_members)
        if family == "nsLTP":
            retained_rso.update(named_ambiguous_nsltp)
        display_tree, kept = main_display_tree(full_tree, retained_rso)
        draw_tree(
            display_tree,
            args.output_dir / f"Figure2{panel}_{family}_phylogeny",
            f"({panel}) {LABELS[family]}: supported-context display",
            primary_members,
            current_rso,
            full=False,
        )
        draw_tree(
            full_tree,
            args.output_dir / f"FigureS1{panel}_{family}_full_phylogeny",
            f"Figure S1{panel}. {LABELS[family]}: complete full-precursor gappyout topology",
            primary_members,
            current_rso,
            full=True,
        )
        manifest_rows.append(
            {
                "panel": panel,
                "family": family,
                "strategy": args.strategy,
                "tree_file": str(tree_path),
                "total_tips": full_tree.count_terminals(),
                "main_display_tips": display_tree.count_terminals(),
                "Rso_primary_members": len(primary_members),
                "display_note": "Main panel retains primary members and named ambiguous nsLTP candidates plus nearest non-RSO context and collapses branches failing the joint threshold; all audited boundary controls and the full uncollapsed topology are in Figure S1",
            }
        )
    with (args.output_dir / "Figure2_phylogeny_panel_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print("PASS representative phylogeny panels")


if __name__ == "__main__":
    main()
