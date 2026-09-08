"""Key-branch circular figures from frozen trees, with separate support audit panels."""

import argparse
import copy
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

from plot_circular_phylogenies import (
    FAMILIES, TITLES, aliases_for, fingerprint, layout, read_tsv,
    save_figure, sha256, stable_groups, write_tsv, xy,
)
from export_phylogeny_itol_package import make_main_display_tree, support_is_joint, tip_category


def map_groups(tree, groups, branch_rows, family):
    records = []
    all_ids = {t.name for t in tree.get_terminals()}
    for group in groups:
        wanted = {"RSO_" + x for x in group["members"].split(";")}
        stable = group["identical_all_taxon_split_strategies"] == 4
        stable_splits = {}
        for row in branch_rows:
            if row["family"] == family and row["rso_group"] == group["members"] and row["joint_support"] == "TRUE":
                stable_splits.setdefault(row["canonical_unrooted_split"], set()).add(row["strategy"])
        if stable:
            candidates = [key for key, strategies in stable_splits.items() if len(strategies) == 4]
            assert candidates
            target_ids = set(min(candidates, key=lambda key: (len(key.split(";")), key)).split(";"))
        else:
            target_ids = wanted
        node = tree.common_ancestor([t for t in tree.get_terminals() if t.name in target_ids])
        ids = {t.name for t in node.get_terminals()}
        if stable:
            assert ids == target_ids
        assert {x for x in ids if x.startswith("RSO_")} == wanted
        canonical = min((tuple(sorted(ids)), tuple(sorted(all_ids - ids))), key=lambda x: (len(x), x))
        strategies = {r["strategy"] for r in branch_rows if r["family"] == family
                      and r["canonical_unrooted_split"] == ";".join(canonical)
                      and r["joint_support"] == "TRUE"}
        if stable:
            assert len(strategies) == 4 and support_is_joint(node.name)
        records.append({**group, "source_support": node.name, "total_sequences": len(ids),
                        "Rso_sequences": len(wanted), "all_taxon_ids": ";".join(sorted(ids)),
                        "source_split_supported_strategies": len(strategies),
                        "highlight_in_main": stable})
    return records


def main_display(tree, primary, ambiguous, records):
    baseline = make_main_display_tree(tree, primary | ambiguous)
    keep = {t.name for t in baseline.get_terminals()}
    for record in records:
        if record["highlight_in_main"]:
            keep.update(record["all_taxon_ids"].split(";"))
    result = copy.deepcopy(tree)
    for tip in list(result.get_terminals()):
        if tip.name not in keep:
            result.prune(tip)
    for node in list(result.get_nonterminals(order="postorder")):
        if node is not result.root and not support_is_joint(node.name):
            result.collapse(node)
    result.ladderize()
    for record in records:
        if record["highlight_in_main"]:
            ids = set(record["all_taxon_ids"].split(";"))
            node = result.common_ancestor([t for t in result.get_terminals() if t.name in ids])
            assert {t.name for t in node.get_terminals()} == ids
            assert node.name == record["source_support"]
    return result


def draw(ax, tree, aliases, primary, records, title, font_size, full=False, audit=False):
    leaves, angles, radii = layout(tree)
    parent_of = {c: p for p in tree.get_nonterminals() for c in p.clades}
    colors, key_nodes, unstable_nodes = {}, [], []
    tip_ids = {t.name for t in leaves}
    for record in records:
        ids = set(record["all_taxon_ids"].split(";"))
        if not ids <= tip_ids:
            continue
        node = tree.common_ancestor([t for t in leaves if t.name in ids])
        if {t.name for t in node.get_terminals()} != ids:
            continue
        if record["highlight_in_main"]:
            assert node.name == record["source_support"]
            for descendant in node.find_clades():
                colors[descendant] = record["color"]
            key_nodes.append((node, record))
        elif full:
            unstable_nodes.append((node, record))

    for parent in tree.get_nonterminals():
        child_angles = [angles[c] for c in parent.clades]
        arc = np.linspace(min(child_angles), max(child_angles), 100)
        ax.plot(*xy(radii[parent], arc), color=colors.get(parent, "#B0B3B5"), lw=.65)
        for child in parent.clades:
            weak = not child.is_terminal() and not support_is_joint(child.name)
            ax.plot(*xy(np.array([radii[parent], radii[child]]), angles[child]),
                    color=colors.get(child, "#B0B3B5"), lw=.8 if child in colors else .55,
                    linestyle=(0, (2, 2)) if full and weak else "-")

    tip_artists = []
    step = math.radians(354 / len(leaves))
    for tip in leaves:
        category = tip_category(tip.name, primary)
        angle = angles[tip]
        degrees = math.degrees(angle) % 360
        left = 90 < degrees < 270
        label = aliases[tip.name]
        if category == "excluded_Rso":
            label = "x:" + label
        artist = ax.text(*xy(1.035, angle), label, ha="right" if left else "left", va="center",
                         rotation=degrees + 180 if left else degrees, rotation_mode="anchor",
                         fontsize=font_size, fontweight="bold" if category == "primary_Rso" else "normal",
                         color=colors.get(tip, "#333333" if category == "primary_Rso" else "#777777"))
        tip_artists.append(artist)

    for node, record in key_nodes:
        a = [angles[t] for t in node.get_terminals()]
        arc = np.linspace(min(a) - step * .4, max(a) + step * .4, 80)
        ax.plot(*xy(1.005, arc), color=record["color"], lw=2.1)

    ax.set(xlim=(-1.55, 1.55), ylim=(-1.55, 1.55), aspect="equal")
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=font_size + 2, fontweight="bold", pad=0)
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    occupied = [a.get_window_extent(renderer).expanded(1.03, 1.08) for a in tip_artists]
    number_artists = []
    if audit:
        for node in tree.get_nonterminals():
            if node not in parent_of or not node.name or "/" not in node.name:
                continue
            radius = (radii[parent_of[node]] + radii[node]) / 2
            degrees = math.degrees(angles[node]) % 360
            rotation = degrees if degrees <= 90 or degrees >= 270 else degrees + 180
            artist = ax.text(*xy(radius, angles[node]), node.name, ha="center", va="center",
                             rotation=rotation, rotation_mode="anchor", fontsize=5.3,
                             color="#555555", bbox=dict(facecolor="white", edgecolor="none", pad=.15))
            number_artists.append(artist)
        occupied.extend(a.get_window_extent(renderer).expanded(1.05, 1.1) for a in number_artists)

    annotations = []
    for node, record in key_nodes + unstable_nodes:
        stable = record["highlight_in_main"]
        color = record["color"] if stable else "#777777"
        radius = (radii[parent_of[node]] + radii[node]) / 2
        anchor = xy(radius, angles[node])
        name = f"{record['group_label']} (n={record['total_sequences']})" if stable else "L5* (unstable)"
        label = name if audit else name + "\n" + record["source_support"]
        artist = ax.text(*anchor, label, ha="center", va="center", fontsize=font_size * .91,
                         fontweight="bold" if stable else "normal", color=color, linespacing=1.25,
                         bbox=dict(facecolor="white", edgecolor="none", pad=.65), zorder=8)
        placed = False
        for dr in (0, -.07, -.14, -.21, -.28, -.35, -.42, -.49, .05):
            for da in (0, 5, -5, 10, -10, 18, -18, 28, -28, 40, -40):
                artist.set_position(xy(max(.1, radius + dr), angles[node] + math.radians(da)))
                box = artist.get_window_extent(renderer).expanded(1.15, 1.18)
                if not any(box.overlaps(b) for b in occupied):
                    occupied.append(box)
                    placed = True
                    break
            if placed:
                break
        if not placed:
            raise ValueError(f"Cannot place {title} {name}")
        ax.scatter(*anchor, s=9 if full else 12, color=color, zorder=7)
        ax.annotate("", xy=anchor, xytext=artist.get_position(),
                    arrowprops=dict(arrowstyle="-", color=color, linewidth=.55), zorder=6)
        annotations.append(artist)
    artists = tip_artists + annotations + number_artists
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    for artist in artists:
        if not all(ax.figure.bbox.contains(*p) for p in artist.get_window_extent(renderer).get_points()):
            raise ValueError(f"Clipped label: {title} {artist.get_text()}")
    return {"tips": len(leaves), "key_groups": [r["group_label"] for n, r in key_nodes],
            "numeric_labels": len(number_artists) if audit else len(annotations),
            "weak_group_reported": bool(unstable_nodes), "clipped_labels": 0}


def footer(fig, size):
    fig.text(.5, .045, "Bold: R. soongarica catalogue   x: audited exclusion   n: sequences in the highlighted group",
             ha="center", fontsize=size, color="#555555")
    fig.text(.5, .02, "SH-aLRT / UFBoot (%)   |   Circular cladogram; radial lengths not to scale; center is not a biological root",
             ha="center", fontsize=size, color="#555555")


def support_audit(tree, aliases, primary, records, title):
    leaves = tree.get_terminals()
    y = {tip: index for index, tip in enumerate(leaves)}
    heights = {}
    for node in tree.find_clades(order="postorder"):
        if node.is_terminal():
            heights[node] = 0
        else:
            heights[node] = 1 + max(heights[c] for c in node.clades)
            values = [y[t] for t in node.get_terminals()]
            y[node] = (min(values) + max(values)) / 2
    maximum = heights[tree.root]
    x = {node: maximum - height for node, height in heights.items()}
    colors = {}
    for record in records:
        if record["highlight_in_main"]:
            ids = set(record["all_taxon_ids"].split(";"))
            node = tree.common_ancestor([t for t in leaves if t.name in ids])
            colors.update({c: record["color"] for c in node.find_clades()})
    fig = plt.figure(figsize=((maximum + 7) * .46, len(leaves) * .16 + 2))
    ax = fig.add_axes([.035, .04, .93, .92])
    numbers = []
    for parent in tree.get_nonterminals():
        ys = [y[c] for c in parent.clades]
        ax.plot([x[parent], x[parent]], [min(ys), max(ys)], color=colors.get(parent, "#999999"), lw=.55)
        for child in parent.clades:
            weak = not child.is_terminal() and not support_is_joint(child.name)
            ax.plot([x[parent], x[child]], [y[child], y[child]], color=colors.get(child, "#999999"),
                    lw=.6, linestyle=(0, (2, 2)) if weak else "-")
            if not child.is_terminal() and child.name and "/" in child.name:
                numbers.append(ax.text(x[child] - .07, y[child] - .10, child.name, fontsize=6.6,
                                       ha="right", va="bottom", color="#333333",
                                       bbox=dict(facecolor="white", edgecolor="none", pad=.15)))
    for tip in leaves:
        category = tip_category(tip.name, primary)
        label = ("x:" if category == "excluded_Rso" else "") + aliases[tip.name]
        ax.text(x[tip] + .2, y[tip], label, ha="left", va="center", fontsize=8.5,
                color=colors.get(tip, "#444444"), fontweight="bold" if category == "primary_Rso" else "normal")
    ax.set(xlim=(-.5, maximum + 6), ylim=(len(leaves), -2))
    ax.axis("off")
    ax.set_title(title + " | all available SH-aLRT / UFBoot (%)", loc="left", fontsize=12, fontweight="bold")
    fig.text(.5, .014, "Support audit only. Cladogram lengths are not evolutionary distances; the display origin is not a biological root.",
             ha="center", fontsize=8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [a.get_window_extent(renderer).expanded(1.02, 1.02) for a in numbers]
    assert not any(a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1:]), "Overlapping audit support labels"
    return fig, {"tips": len(leaves), "numeric_labels": len(numbers), "layout": "rectangular support audit",
                 "key_groups": [r["group_label"] for r in records if r["highlight_in_main"]], "clipped_labels": 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output, revision = args.output, args.revision
    data, qa = output / "04_Newick_and_annotations", output / "05_QA"
    data.mkdir(parents=True, exist_ok=True)
    qa.mkdir(parents=True, exist_ok=True)
    catalogue = read_tsv(revision / "03_annotation_rescue/final_amp_primary_catalogue.tsv")
    primary = {f: {r["member_id"] for r in catalogue if r["family"] == f} for f in FAMILIES}
    ambiguous = {r["member_id"] for r in read_tsv(revision / "01_nslTP_curation/reclassified_candidates.tsv")
                 if r["member_id"] and r["final_classification"] == "B_AMBIGUOUS_NSLTP_PROLAMIN_BOUNDARY"}
    groups = stable_groups(revision, primary)
    branch_rows = read_tsv(revision / "02_phylogeny/phylogeny_branch_support.tsv")
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "svg.fonttype": "none"})
    proofs, all_records, label_rows, panels = [], [], [], []
    for panel, family in zip("ABC", FAMILIES):
        source = revision / f"02_phylogeny/{family}/{family}.full.gappyout.treefile"
        tree = Phylo.read(source, "newick")
        before = fingerprint(tree)
        tree.ladderize()
        records = map_groups(tree, [g for g in groups if g["family"] == family], branch_rows, family)
        display = main_display(tree, primary[family], ambiguous if family == "nsLTP" else set(), records)
        aliases = aliases_for(tree, primary[family], family)
        all_records.extend(records)
        for full, audit, t, name in (
            (False, False, display, f"Figure2{panel}_{family}_key_clades"),
            (True, False, tree, f"FigureS1{panel}_{family}_key_clades"),
            (True, True, tree, f"Support_audit_{panel}_{family}_all_nodes"),
        ):
            title = f"({panel}) {TITLES[family]} | {len(t.get_terminals())} sequences"
            if audit:
                fig, proof = support_audit(t, aliases, primary[family], records, title)
            else:
                fig = plt.figure(figsize=(13, 13) if full else (9, 9))
                ax = fig.add_axes([.02, .095, .96, .85])
                proof = draw(ax, t, aliases, primary[family], records, title, 8.4 if full else 11, full)
                footer(fig, 8.4 if full else 7.5)
            save_figure(fig, name, output / ("03_Figures/Support_audit" if audit else "03_Figures"))
            plt.close(fig)
            proofs.append({"file": name, **proof, "full_tree": full, "audit_all_support": audit})
            if not audit:
                Phylo.write(t, data / f"{name}.nwk", "newick", format_branch_length="%1.10f")
                (data / f"{name}.itol_labels.txt").write_text(
                    "LABELS\nSEPARATOR TAB\nDATA\n" + "\n".join(
                        f"{tip.name}\t{'x:' if tip_category(tip.name, primary[family]) == 'excluded_Rso' else ''}{aliases[tip.name]}"
                        for tip in t.get_terminals()) + "\n", encoding="utf-8")
        for tip in tree.get_terminals():
            label_rows.append({"family": family, "newick_id": tip.name, "display_label": aliases[tip.name],
                               "category": tip_category(tip.name, primary[family]),
                               "shown_in_main": tip.name in {t.name for t in display.get_terminals()}})
        assert fingerprint(tree) == before
        (data / source.name.replace(".treefile", ".original.treefile")).write_bytes(source.read_bytes())
        panels.append((display, aliases, primary[family], records, panel, family))
    fig = plt.figure(figsize=(7, 9.3))
    for position, (tree, aliases, members, records, panel, family) in zip(
        ((.025, .57, .46, .39), (.515, .57, .46, .39), (.10, .105, .80, .46)), panels
    ):
        ax = fig.add_axes(position)
        proof = draw(ax, tree, aliases, members, records, f"({panel}) {TITLES[family]}",
                     7.2 if panel != "C" else 7.7)
        proofs.append({"file": "Figure2_combined_key_clades", "panel": panel, **proof})
    fig.text(.5, .045, "Bold: R. soongarica catalogue   x: audited exclusion", ha="center", fontsize=7.2)
    fig.text(.5, .022, "n: sequences in group   |   SH-aLRT / UFBoot (%)", ha="center", fontsize=7.2)
    save_figure(fig, "Figure2_combined_key_clades", output / "03_Figures")
    plt.close(fig)
    assert sum(r["highlight_in_main"] for r in all_records) == 7
    write_tsv(data / "key_branch_evidence.tsv", all_records)
    write_tsv(data / "tip_label_key.tsv", label_rows)
    originals = data / "All_12_original_ML_trees"
    originals.mkdir(exist_ok=True)
    for family in FAMILIES:
        for region in ("full", "core"):
            for trim in ("untrimmed", "gappyout"):
                path = revision / f"02_phylogeny/{family}/{family}.{region}.{trim}.treefile"
                (originals / path.name).write_bytes(path.read_bytes())
    (data / "all_numeric_branch_support.tsv").write_bytes((revision / "02_phylogeny/phylogeny_branch_support.tsv").read_bytes())
    (qa / "key_figure_checks.json").write_text(json.dumps(proofs, indent=2), encoding="utf-8")
    print(json.dumps(proofs))


if __name__ == "__main__":
    main()
