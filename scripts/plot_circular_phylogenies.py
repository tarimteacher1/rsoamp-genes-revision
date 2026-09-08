"""Circular review figures from frozen ML trees; no inference or rerooting."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from matplotlib.lines import Line2D

from export_phylogeny_itol_package import make_main_display_tree, support_is_joint, tip_category


FAMILIES = ("Defensin", "Snakin_GASA", "nsLTP")
TITLES = {"Defensin": "Defensin", "Snakin_GASA": "Snakin/GASA", "nsLTP": "nsLTP"}
ORIGINS = {
    "primary_Rso": ("R. soongarica", "#B64043"),
    "excluded_Rso": ("Rso exclusion", "#989898"),
    "Tamarix": ("T. austromongolica", "#237B8D"),
    "positive_reference": ("Reviewed reference", "#444444"),
    "boundary_control": ("2S/prolamin control", "#AB7A22"),
}
GROUP_COLORS = ("#287CB4", "#CE6941", "#26947B", "#BA4D8E", "#798832")


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(tree):
    return sorted(
        (tuple(sorted(t.name for t in c.get_terminals())), c.name or "", c.branch_length or 0.0)
        for c in tree.find_clades()
    )


def stable_groups(revision, primary):
    groups = read_tsv(revision / "02_phylogeny/rso_supported_group_sensitivity.tsv")
    branches = read_tsv(revision / "02_phylogeny/phylogeny_branch_support.tsv")
    selected = []
    for family, prefix in zip(FAMILIES, ("D", "G", "L")):
        rows = []
        for row in groups:
            members = set(row["rso_group"].split(";"))
            if row["family"] != family or len(members) < 2:
                continue
            if not members <= primary[family] or int(row["supported_strategy_count"]) < 3:
                continue
            exact = defaultdict(set)
            for branch in branches:
                if (branch["family"] == family and branch["rso_group"] == row["rso_group"]
                        and branch["joint_support"] == "TRUE"):
                    exact[branch["canonical_unrooted_split"]].add(branch["strategy"])
            full_count = max(map(len, exact.values()), default=0)
            rows.append({
                "family": family,
                "members": row["rso_group"],
                "rso_restricted_strategies": int(row["supported_strategy_count"]),
                "identical_all_taxon_split_strategies": full_count,
                "color": "",
                "group_label": "",
                "interpretation": "all-taxon split stable" if full_count >= 3 else "Rso-restricted only; not all-taxon stable",
            })
        rows.sort(key=lambda x: (-x["identical_all_taxon_split_strategies"], x["members"]))
        for index, row in enumerate(rows):
            row["color"] = GROUP_COLORS[index % len(GROUP_COLORS)]
            row["group_label"] = f"{prefix}{index + 1}" + ("*" if row["identical_all_taxon_split_strategies"] < 3 else "")
        selected.extend(rows)
    return selected


def aliases_for(tree, primary, family):
    aliases = {}
    tamarix_index = 0
    family_short = {"Defensin": "DEF", "Snakin_GASA": "GAS", "nsLTP": "LTP"}[family]
    for name in sorted(t.name for t in tree.get_terminals()):
        if name.startswith("RSO_"):
            label = name[4:]
            if not label.startswith("Rs"):
                label = f"RsoX{sum(v.startswith('RsoX') for v in aliases.values()) + 1:02d}"
        elif name.startswith("TAU_"):
            tamarix_index += 1
            label = f"Ta{family_short}{tamarix_index:02d}"
        elif name.startswith("AMBPRO_"):
            label = "Pro:" + name[7:].removeprefix("AMB").removeprefix("_")
        elif name.startswith("NEG2S_"):
            label = "2S:" + name[6:].replace("NEG2S_", "")
        else:
            label = name.removeprefix("UNI_").removeprefix("NSPOS_").removeprefix("POS_")
        aliases[name] = label
    if len(set(aliases.values())) != len(aliases):
        raise ValueError("Display aliases are not unique")
    assert all(aliases.values())
    return aliases


def layout(tree):
    leaves = tree.get_terminals()
    angles = {tip: math.radians(93 + i * 354 / len(leaves)) for i, tip in enumerate(leaves)}
    heights = {}
    indices = {tip: i for i, tip in enumerate(leaves)}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            heights[node] = 0
        else:
            descendants = node.get_terminals()
            lo, hi = min(indices[t] for t in descendants), max(indices[t] for t in descendants)
            angles[node] = (angles[leaves[lo]] + angles[leaves[hi]]) / 2
            heights[node] = 1 + max(heights[child] for child in node.clades)
    maximum = heights[tree.root]
    radii = {node: 0.06 + 0.92 * (1 - height / maximum) for node, height in heights.items()}
    return leaves, angles, radii


def xy(radius, angle):
    return radius * np.cos(angle), radius * np.sin(angle)


def draw_circle(ax, tree, aliases, primary, groups, title, full, font_size):
    leaves, angles, radii = layout(tree)
    branch_group = {}
    group_nodes = []
    support_artists = []
    group_artists = []
    for group in groups:
        wanted = {"RSO_" + x for x in group["members"].split(";")}
        found = [tip for tip in leaves if tip.name in wanted]
        if len(found) != len(wanted):
            continue
        ancestor = tree.common_ancestor(found)
        present = {tip.name[4:] for tip in ancestor.get_terminals() if tip.name.startswith("RSO_")}
        if present != set(group["members"].split(";")):
            continue
        for node in ancestor.find_clades():
            branch_group[node] = group
        group_nodes.append((ancestor, group))

    for parent in tree.get_nonterminals():
        child_angles = [angles[c] for c in parent.clades]
        arc = np.linspace(min(child_angles), max(child_angles), max(8, int(math.degrees(max(child_angles) - min(child_angles))) * 2))
        group = branch_group.get(parent)
        color = group["color"] if group else "#777777"
        ax.plot(*xy(radii[parent], arc), color=color, lw=0.65, solid_capstyle="round")
        for child in parent.clades:
            group = branch_group.get(child)
            color = group["color"] if group else "#777777"
            weak = not child.is_terminal() and not support_is_joint(child.name)
            if full and weak:
                color = "#B8B8B8"
            ax.plot(*xy(np.array([radii[parent], radii[child]]), angles[child]),
                    color=color, lw=0.68 if full else 0.95,
                    linestyle=(0, (2, 1.7)) if full and weak else "-", solid_capstyle="round")
            if full and not child.is_terminal() and support_is_joint(child.name):
                radius = (radii[parent] + radii[child]) / 2
                degrees = math.degrees(angles[child]) % 360
                rotation = degrees if degrees <= 90 or degrees >= 270 else degrees + 180
                artist = ax.text(*xy(radius, angles[child]), child.name, ha="center", va="center",
                        rotation=rotation, rotation_mode="anchor", fontsize=max(4.8, font_size * 0.61),
                        color="#333333", bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.2, "alpha": 0.9})
                support_artists.append(artist)

    label_artists = []
    step = math.radians(354 / len(leaves))
    for tip in leaves:
        angle = angles[tip]
        category = tip_category(tip.name, primary)
        color = ORIGINS[category][1]
        arc = np.linspace(angle - step * 0.37, angle + step * 0.37, 10)
        ax.plot(*xy(1.015, arc), color=color, lw=3.0 if full else 3.5, solid_capstyle="butt")
        degrees = math.degrees(angle) % 360
        left = 90 < degrees < 270
        artist = ax.text(*xy(1.055, angle), aliases[tip.name], ha="right" if left else "left", va="center",
                        rotation=degrees + 180 if left else degrees, rotation_mode="anchor",
                        fontsize=font_size, fontweight="bold" if category == "primary_Rso" else "normal",
                        color="#262626" if category != "excluded_Rso" else "#777777")
        label_artists.append(artist)

    for ancestor, group in group_nodes:
        a = [angles[t] for t in ancestor.get_terminals()]
        arc = np.linspace(min(a) - step * .4, max(a) + step * .4, 60)
        radius = 0.994
        tentative = group["identical_all_taxon_split_strategies"] < 3
        ax.plot(*xy(radius, arc), color=group["color"], lw=2.2,
                linestyle=(0, (2, 1.5)) if tentative else "-")
        label_radius = max(.18, radii[ancestor] - .075)
        artist = ax.text(*xy(label_radius, angles[ancestor]), group["group_label"],
                color=group["color"], fontsize=font_size * .85, fontweight="bold",
                ha="center", va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": .35})
        group_artists.append((artist, label_radius, angles[ancestor]))

    ax.set(xlim=(-1.63, 1.63), ylim=(-1.63, 1.63), aspect="equal")
    ax.axis("off")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=font_size + 2, pad=0)
    if full:
        # Separate group names from numeric support without moving any tree nodes.
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
        occupied = [a.get_window_extent(renderer).expanded(1.1, 1.2) for a in support_artists]
        for artist, radius, angle in group_artists:
            placed = False
            for radial_offset in (0, .045, .09, .135, .18):
                for angle_offset in (0, 4, -4, 8, -8):
                    artist.set_position(xy(max(.12, radius - radial_offset), angle + math.radians(angle_offset)))
                    box = artist.get_window_extent(renderer).expanded(1.25, 1.3)
                    if not any(box.overlaps(other) for other in occupied):
                        occupied.append(box)
                        placed = True
                        break
                if placed:
                    break
            if not placed:
                raise ValueError(f"Cannot separate group/support labels: {artist.get_text()}")
    return label_artists


def legend(fig, size=8, y=.02):
    handles = [Line2D([0], [0], color=color, lw=4, label=label) for label, color in ORIGINS.values()]
    fig.legend(handles=handles, ncol=3, loc="lower center", bbox_to_anchor=(.5, y),
               frameon=False, fontsize=size, handlelength=1.1, columnspacing=1.2)


def save_figure(fig, base, root):
    for extension, folder in (("pdf", "Vector_PDF"), ("svg", "Editable_SVG"), ("png", "High_Resolution_PNG")):
        output = root / folder / f"{base}.{extension}"
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=600, facecolor="white")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fig_root = args.output / "03_Figures"
    data_root = args.output / "04_Newick_and_annotations"
    data_root.mkdir(parents=True, exist_ok=True)
    qa_root = args.output / "05_QA"
    qa_root.mkdir(parents=True, exist_ok=True)
    catalogue = read_tsv(args.revision / "03_annotation_rescue/final_amp_primary_catalogue.tsv")
    primary = {family: {r["member_id"] for r in catalogue if r["family"] == family} for family in FAMILIES}
    ambiguous = {r["member_id"] for r in read_tsv(args.revision / "01_nslTP_curation/reclassified_candidates.tsv")
                 if r["member_id"] and r["final_classification"] == "B_AMBIGUOUS_NSLTP_PROLAMIN_BOUNDARY"}
    groups = stable_groups(args.revision, primary)
    write_tsv(data_root / "colored_group_definitions.tsv", groups)
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9})
    manifest, labels, proofs, panels = [], [], [], []
    for panel, family in zip("ABC", FAMILIES):
        source = args.revision / f"02_phylogeny/{family}/{family}.full.gappyout.treefile"
        full_tree = Phylo.read(source, "newick")
        original = fingerprint(full_tree)
        full_tree.ladderize()
        main_tree = make_main_display_tree(full_tree, primary[family] | (ambiguous if family == "nsLTP" else set()))
        aliases = aliases_for(full_tree, primary[family], family)
        family_groups = [g for g in groups if g["family"] == family]
        for tip in full_tree.get_terminals():
            labels.append({"family": family, "newick_id": tip.name, "display_label": aliases[tip.name],
                           "category": tip_category(tip.name, primary[family]), "shown_in_main": tip.name in {t.name for t in main_tree.get_terminals()}})
        for full, tree, base in ((False, main_tree, f"Figure2{panel}_{family}_circular"),
                                 (True, full_tree, f"FigureS1{panel}_{family}_circular_full")):
            fig = plt.figure(figsize=(13, 13) if full else (9, 9))
            ax = fig.add_axes([.015, .09, .97, .87])
            artists = draw_circle(ax, tree, aliases, primary[family], family_groups,
                                  f"({panel}) {TITLES[family]}  |  {tree.count_terminals()} sequences", full,
                                  font_size=8.4 if full else 11.5)
            legend(fig, size=9 if full else 8.5, y=.025)
            fig.text(.5, .012, "Circular cladogram; branch lengths not to scale. Center does not imply an evolutionary root.",
                     fontsize=8.5 if full else 7.5, ha="center")
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            bounds = fig.bbox
            clipped = [a.get_text() for a in artists if not bounds.contains(*a.get_window_extent(renderer).get_points()[0])
                       or not bounds.contains(*a.get_window_extent(renderer).get_points()[1])]
            if clipped:
                raise ValueError(f"Clipped labels in {base}: {clipped}")
            save_figure(fig, base, fig_root)
            plt.close(fig)
            Phylo.write(tree, data_root / f"{base}.nwk", "newick", format_branch_length="%1.10f")
            tip_names = [tip.name for tip in tree.get_terminals()]
            (data_root / f"{base}.itol_labels.txt").write_text(
                "LABELS\nSEPARATOR TAB\nDATA\n" + "\n".join(f"{name}\t{aliases[name]}" for name in tip_names) + "\n",
                encoding="utf-8")
            (data_root / f"{base}.itol_origin_strip.txt").write_text(
                "DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tSequence origin\nCOLOR\t#444444\nSTRIP_WIDTH\t18\nDATA\n"
                + "\n".join(f"{name}\t{ORIGINS[tip_category(name, primary[family])][1]}\t{ORIGINS[tip_category(name, primary[family])][0]}"
                             for name in tip_names) + "\n", encoding="utf-8")
            proofs.append({"figure": base, "tips": tree.count_terminals(), "blank_tip_labels": 0, "clipped_tip_labels": clipped,
                           "layout": "circular cladogram", "rerooted": False, "main_weak_nodes_collapsed": not full})
        assert fingerprint(full_tree) == original, "Plotting modified input topology, support, or branch lengths"
        source_copy = data_root / f"{family}.full.gappyout.original.treefile"
        source_copy.write_bytes(source.read_bytes())
        manifest.append({"family": family, "source": str(source.resolve()), "source_sha256": sha256(source),
                         "full_tips": full_tree.count_terminals(), "main_tips": main_tree.count_terminals(),
                         "original_tree_copy_sha256": sha256(source_copy), "topology_support_lengths_preserved": True})
        panels.append((main_tree, aliases, primary[family], family_groups, panel, family))

    combined = plt.figure(figsize=(7, 9.3))
    positions = ((.025, .57, .46, .39), (.515, .57, .46, .39), (.14, .11, .72, .45))
    combined_labels = []
    for position, (tree, aliases, members, family_groups, panel, family) in zip(positions, panels):
        ax = combined.add_axes(position)
        combined_labels.extend(draw_circle(ax, tree, aliases, members, family_groups,
                    f"({panel}) {TITLES[family]}", False, 7.2 if panel != "C" else 7.7))
        ax.set(xlim=(-1.72, 1.72), ylim=(-1.72, 1.72))
    legend(combined, size=7.2, y=.018)
    combined.canvas.draw()
    renderer = combined.canvas.get_renderer()
    for artist in combined_labels:
        if not all(combined.bbox.contains(*point) for point in artist.get_window_extent(renderer).get_points()):
            raise ValueError(f"Clipped combined-panel label: {artist.get_text()}")
    save_figure(combined, "Figure2_combined_circular", fig_root)
    plt.close(combined)
    write_tsv(data_root / "tip_label_key.tsv", labels)
    write_tsv(data_root / "figure_source_manifest.tsv", manifest)
    all_trees = data_root / "All_12_original_ML_trees"
    all_trees.mkdir(exist_ok=True)
    for family in FAMILIES:
        for region in ("full", "core"):
            for trim in ("untrimmed", "gappyout"):
                original_path = args.revision / f"02_phylogeny/{family}/{family}.{region}.{trim}.treefile"
                (all_trees / original_path.name).write_bytes(original_path.read_bytes())
    (data_root / "all_numeric_branch_support.tsv").write_bytes((args.revision / "02_phylogeny/phylogeny_branch_support.tsv").read_bytes())
    (qa_root / "circular_figure_checks.json").write_text(json.dumps(proofs, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"figures": 7, "main_tips": [x["main_tips"] for x in manifest], "full_tips": [x["full_tips"] for x in manifest],
                      "blank_labels": 0, "source_topology_support_lengths_preserved": True}))


if __name__ == "__main__":
    main()
