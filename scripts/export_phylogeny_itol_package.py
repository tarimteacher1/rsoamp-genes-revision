#!/usr/bin/env python3
"""Export reproducible Newick trees and iTOL annotations for the revision figures."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import re
import shutil
from pathlib import Path

from Bio import Phylo


FAMILIES = ["Defensin", "Snakin_GASA", "nsLTP"]
PANELS = {"Defensin": "A", "Snakin_GASA": "B", "nsLTP": "C"}
DISPLAY_NAMES = {
    "Defensin": "Defensin",
    "Snakin_GASA": "Snakin_GASA",
    "nsLTP": "nsLTP",
}
STRATEGIES = ["full.untrimmed", "full.gappyout", "core.untrimmed", "core.gappyout"]
COLORS = {
    "primary_Rso": "#B33A3A",
    "excluded_Rso": "#D47A2C",
    "Tamarix": "#2E6F9E",
    "positive_reference": "#555555",
    "boundary_control": "#8C6D9E",
}
CATEGORY_LABELS = {
    "primary_Rso": "R. soongarica primary catalogue",
    "excluded_Rso": "R. soongarica audited exclusion",
    "Tamarix": "Tamarix austromongolica",
    "positive_reference": "Reviewed family reference",
    "boundary_control": "Prolamin/2S boundary control",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def support_is_joint(text: str | None) -> bool:
    if not text:
        return False
    match = re.fullmatch(r"([0-9.]+)/([0-9.]+)", str(text))
    if not match:
        return False
    sh_alrt, ufboot = map(float, match.groups())
    return sh_alrt >= 80 and ufboot >= 95


def tip_category(name: str, primary_members: set[str]) -> str:
    if name.startswith("RSO_"):
        return "primary_Rso" if name.removeprefix("RSO_") in primary_members else "excluded_Rso"
    if name.startswith("TAU_"):
        return "Tamarix"
    if name.startswith(("NEG2S_", "AMBPRO_")):
        return "boundary_control"
    return "positive_reference"


def itol_display_label(name: str, primary_members: set[str]) -> str:
    if name.startswith("RSO_"):
        member = name.removeprefix("RSO_")
        return member if member in primary_members else f"{member}__excluded"
    if name.startswith("TAU_"):
        return name.removeprefix("TAU_")
    if name.startswith("UNI_"):
        return name.removeprefix("UNI_")
    if name.startswith("NSPOS_"):
        return name.removeprefix("NSPOS_")
    if name.startswith("NEG2S_"):
        return name.removeprefix("NEG2S_") + "__2S_control"
    return name


def make_main_display_tree(tree, retained_rso_members: set[str]):
    rso_tips = [
        tip
        for tip in tree.get_terminals()
        if (tip.name or "").startswith("RSO_")
        and (tip.name or "").removeprefix("RSO_") in retained_rso_members
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
        if clade is not root and not support_is_joint(clade.name):
            display.collapse(clade)
    display.ladderize()
    return display


def write_tree(tree, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, path, "newick", format_branch_length="%1.10f")


def write_colorstrip(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tTaxon category",
        "COLOR\t#333333",
        "STRIP_WIDTH\t30",
        "MARGIN\t5",
        "BORDER_WIDTH\t0",
        "SHOW_INTERNAL\t0",
        "LEGEND_TITLE\tTaxon category",
        "LEGEND_SHAPES\t" + "\t".join("1" for _ in CATEGORY_LABELS),
        "LEGEND_COLORS\t" + "\t".join(COLORS[key] for key in CATEGORY_LABELS),
        "LEGEND_LABELS\t" + "\t".join(CATEGORY_LABELS[key] for key in CATEGORY_LABELS),
        "DATA",
    ]
    lines.extend(f"{row['itol_label']}\t{row['color']}\t{row['category_label']}" for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_label_colors(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["TREE_COLORS", "SEPARATOR TAB", "DATA"]
    lines.extend(f"{row['itol_label']}\tlabel\t{row['color']}\tbold" for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def metadata_rows(tree, family: str, primary_members: set[str], manifest_by_id: dict[str, dict[str, str]], simplify: bool):
    rows = []
    seen = set()
    for tip in tree.get_terminals():
        raw_id = tip.name or ""
        label = itol_display_label(raw_id, primary_members) if simplify else raw_id
        if label in seen:
            raise ValueError(f"Non-unique iTOL label after simplification: {family} {label}")
        seen.add(label)
        category = tip_category(raw_id, primary_members)
        source = manifest_by_id.get(raw_id, {})
        rows.append(
            {
                "family": family,
                "raw_taxon_id": raw_id,
                "itol_label": label,
                "category": category,
                "category_label": CATEGORY_LABELS[category],
                "color": COLORS[category],
                "source_group": source.get("source_group", ""),
                "source_id": source.get("source_id", ""),
                "length_aa": source.get("length_aa", ""),
                "note": source.get("note", ""),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256_manifest(output_dir: Path) -> None:
    rows = []
    for path in sorted(p for p in output_dir.rglob("*") if p.is_file() and p.name != "SHA256SUMS.tsv"):
        rows.append(
            {
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "relative_path": path.relative_to(output_dir).as_posix(),
            }
        )
    write_tsv(output_dir / "SHA256SUMS.tsv", rows)


def write_readme(output_dir: Path, manifest_rows: list[dict[str, str]]) -> None:
    mapping = "\n".join(
        f"- Figure 2{row['panel']}: `01_Main_Figure2/{row['main_tree']}`; Figure S1{row['panel']}: `02_Full_FigureS1/{row['full_tree']}`."
        for row in manifest_rows
    )
    readme = f"""# R. soongarica AMP phylogeny iTOL package

## Figure mapping

{mapping}

The three current main and supplementary panels use the `full.gappyout` strategy. The 12 ML trees under `03_All_Sensitivity_ML_Trees` cover full precursor versus mature/core sequence and untrimmed versus trimAl `gappyout` alignment for each family. The matching `.contree` files under `03b_Optional_Consensus_Trees` are included for audit only; they were not the source of the current figures.

## Critical interpretation rules

1. IQ-TREE inferred these as unrooted maximum-likelihood trees. In iTOL, use an unrooted or circular display unless a biologically defensible outgroup is specified. The arbitrary left-hand origin of a rectangular rendering is not an ancestral root.
2. Internal labels in the raw trees are `SH-aLRT/UFBoot`, based on 1,000 replicates for each statistic. Joint support means SH-aLRT >= 80 and UFBoot >= 95.
3. The main-display trees retain primary R. soongarica members, the four named ambiguous nsLTP boundary candidates, and each retained R. soongarica tip's nearest non-RSO context. Nodes failing the joint support threshold were collapsed. Multifurcations therefore represent unresolved relationships and must not be manually converted into supported bifurcations.
4. The complete Figure S1 trees are uncollapsed. Show support labels there, preferably hiding labels below the joint threshold. Main panels may omit numerical labels if the legend states that unsupported nodes were collapsed.
5. Stable groups are sensitivity-qualified supported clades, not formal subfamilies. The primary interpretation uses agreement across all four strategies, not topology from one selected tree alone.
6. Files ending in `_display_labels.nwk` use shortened unique labels for iTOL layout. The raw full and sensitivity trees retain the exact analysis taxon IDs. TSV metadata files provide the mapping.

## iTOL annotation files

For each panel, upload the matching `*_label_colors.txt` and `*_colorstrip.txt` files after loading the Newick tree. The corresponding `*_tip_metadata.tsv` records the raw ID, displayed label, source category, accession/source ID, sequence length, and curation note.

## Recommended layout

- Figure 2: circular or unrooted display; readable leaf labels; primary R. soongarica genes in dark red; audited exclusions in orange; Tamarix in blue; reviewed references in gray; boundary controls in purple.
- Figure S1: circular display is preferred for the 104-tip Defensin and 153-tip nsLTP trees; retain the scale bar and export as PDF/SVG.
- Do not rotate branches in a way that implies a changed topology. Branch rotation around a node is cosmetic, but rerooting changes the visual interpretation and must be reported.
"""
    (output_dir / "README_iTOL.md").write_text(readme, encoding="utf-8")
    readme_zh = f"""# 红砂 AMP 系统发育树 iTOL 重绘包

## 图与树文件对应关系

{mapping}

当前主图和 Figure S1 均以 `full.gappyout` 作为展示策略；`03_All_Sensitivity_ML_Trees` 中另有每个家族的 4 种分析树，共 12 棵：全长/核心区分别配合未剪切/trimAl gappyout 比对。`03b_Optional_Consensus_Trees` 中另附 12 棵 `.contree` 供审计；它们不是当前图件的来源，不能在不修改方法说明的情况下替换主图树。

## 重绘时必须遵守的解释边界

1. 这些是 IQ-TREE 推断的无根最大似然树。iTOL 建议使用无根图或环形图；矩形图左侧起点不能解释为祖先根。没有可靠外群时不要人为定根。
2. 原始树内部节点标签为 `SH-aLRT/UFBoot`，两项均使用 1,000 次重复。联合支持阈值为 SH-aLRT >= 80 且 UFBoot >= 95。
3. 主图显示树保留红砂主目录成员、4 个具名 nsLTP 边界排除候选，以及每个红砂末端最近的非红砂参考序列。未达到联合阈值的节点已折叠。因此三叉或多叉表示“关系未解析”，不是树文件错误，也不能在 iTOL 中手工改成有支持的二叉关系。
4. Figure S1 是未经折叠的完整树，适合显示数值支持度。主图可隐藏节点数值，但图注必须保留“低支持节点已折叠”的说明。
5. 稳定关系只能称为经过 4 策略敏感性检验的 supported clades，不能恢复为正式 subfamilies。正文结论依据至少 3/4 策略重复出现的关系，而不是单棵展示树。
6. 主图 `_display_labels.nwk` 使用便于排版的短标签；完整树和 12 棵敏感性树保留原始分析 ID。TSV 元数据提供两者映射。

## iTOL 注释

载入 Newick 后，再上传同一面板对应的 `*_label_colors.txt` 和 `*_colorstrip.txt`。颜色含义：红色为红砂主目录，橙色为审核后排除的红砂候选，蓝色为 Tamarix austromongolica，灰色为已审校家族参考，紫色为 prolamin/2S 边界对照。

## 建议输出

- Figure 2：环形或无根布局，优先保证叶标签可读；保留比例尺；导出 PDF/SVG。
- Figure S1：Defensin 104 个末端、nsLTP 153 个末端，优先环形布局并展示支持度。
- 围绕内部节点旋转分支只改变排版，不改变拓扑；重新定根会改变视觉解释，必须有明确依据并在方法或图注说明。
"""
    (output_dir / "README_iTOL_中文.md").write_text(readme_zh, encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True)

    primary = read_tsv(args.revision / "03_annotation_rescue/final_amp_primary_catalogue.tsv")
    primary_by_family = {
        family: {row["member_id"] for row in primary if row["family"] == family}
        for family in FAMILIES
    }
    reclassified = read_tsv(args.revision / "01_nslTP_curation/reclassified_candidates.tsv")
    named_ambiguous = {
        row["member_id"]
        for row in reclassified
        if row.get("member_id")
        and row.get("final_classification") == "B_AMBIGUOUS_NSLTP_PROLAMIN_BOUNDARY"
    }

    manifest_rows = []
    for family in FAMILIES:
        panel = PANELS[family]
        family_label = DISPLAY_NAMES[family]
        taxon_manifest = read_tsv(args.revision / f"02_phylogeny/{family}/{family}.taxon_manifest.tsv")
        manifest_by_id = {row["taxon_id"]: row for row in taxon_manifest}
        primary_members = primary_by_family[family]

        selected_path = args.revision / f"02_phylogeny/{family}/{family}.full.gappyout.treefile"
        full_tree = Phylo.read(selected_path, "newick")
        full_tree.ladderize()
        retained = set(primary_members)
        if family == "nsLTP":
            retained.update(named_ambiguous)
        display_tree = make_main_display_tree(full_tree, retained)
        for tip in display_tree.get_terminals():
            tip.name = itol_display_label(tip.name or "", primary_members)

        main_name = f"Figure2{panel}_{family_label}_supported_context_display_labels.nwk"
        full_name = f"FigureS1{panel}_{family_label}_full_gappyout_raw_labels.nwk"
        write_tree(display_tree, args.output_dir / "01_Main_Figure2" / main_name)
        (args.output_dir / "02_Full_FigureS1").mkdir(parents=True, exist_ok=True)
        shutil.copy2(selected_path, args.output_dir / "02_Full_FigureS1" / full_name)

        raw_tree = Phylo.read(selected_path, "newick")
        # Display-tree labels are already simplified, so map through the selected raw tree.
        raw_display = make_main_display_tree(Phylo.read(selected_path, "newick"), retained)
        main_rows = metadata_rows(raw_display, family, primary_members, manifest_by_id, simplify=True)
        full_rows = metadata_rows(raw_tree, family, primary_members, manifest_by_id, simplify=False)

        prefix = f"Figure2{panel}_{family_label}"
        write_tsv(args.output_dir / "04_iTOL_Annotations" / f"{prefix}_tip_metadata.tsv", main_rows)
        write_colorstrip(args.output_dir / "04_iTOL_Annotations" / f"{prefix}_colorstrip.txt", main_rows)
        write_label_colors(args.output_dir / "04_iTOL_Annotations" / f"{prefix}_label_colors.txt", main_rows)
        sprefix = f"FigureS1{panel}_{family_label}"
        write_tsv(args.output_dir / "04_iTOL_Annotations" / f"{sprefix}_tip_metadata.tsv", full_rows)
        write_colorstrip(args.output_dir / "04_iTOL_Annotations" / f"{sprefix}_colorstrip.txt", full_rows)
        write_label_colors(args.output_dir / "04_iTOL_Annotations" / f"{sprefix}_label_colors.txt", full_rows)

        for strategy in STRATEGIES:
            source = args.revision / f"02_phylogeny/{family}/{family}.{strategy}.treefile"
            target = args.output_dir / "03_All_Sensitivity_ML_Trees" / family / f"{family_label}.{strategy}.raw_labels.nwk"
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            consensus_source = args.revision / f"02_phylogeny/{family}/{family}.{strategy}.contree"
            consensus_target = args.output_dir / "03b_Optional_Consensus_Trees" / family / f"{family_label}.{strategy}.consensus.raw_labels.nwk"
            consensus_target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(consensus_source, consensus_target)

        manifest_rows.append(
            {
                "panel": panel,
                "family": family,
                "selected_strategy": "full.gappyout",
                "main_tree": main_name,
                "main_tips": str(display_tree.count_terminals()),
                "full_tree": full_name,
                "full_tips": str(full_tree.count_terminals()),
                "primary_Rso_members": str(len(primary_members)),
                "named_ambiguous_Rso_members": str(len(named_ambiguous) if family == "nsLTP" else 0),
            }
        )

    write_tsv(args.output_dir / "TREE_MANIFEST.tsv", manifest_rows)
    (args.output_dir / "05_Supporting_Tables").mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.revision / "02_phylogeny/alignment_comparison.tsv", args.output_dir / "05_Supporting_Tables/alignment_comparison.tsv")
    shutil.copy2(args.revision / "02_phylogeny/phylogeny_sensitivity_summary.tsv", args.output_dir / "05_Supporting_Tables/phylogeny_sensitivity_summary.tsv")
    shutil.copy2(args.revision / "02_phylogeny/rso_supported_group_sensitivity.tsv", args.output_dir / "05_Supporting_Tables/rso_supported_group_sensitivity.tsv")
    shutil.copy2(args.revision / "02_phylogeny/phylogeny_branch_support.tsv", args.output_dir / "05_Supporting_Tables/phylogeny_branch_support.tsv")
    shutil.copy2(args.revision / "02_phylogeny/phylogeny_sensitivity_summary.md", args.output_dir / "05_Supporting_Tables/phylogeny_sensitivity_summary.md")
    shutil.copy2(args.revision / "02_phylogeny/software_versions.txt", args.output_dir / "05_Supporting_Tables/software_versions.txt")
    write_readme(args.output_dir, manifest_rows)
    sha256_manifest(args.output_dir)
    print(f"PASS exported {args.output_dir}")


if __name__ == "__main__":
    main()
