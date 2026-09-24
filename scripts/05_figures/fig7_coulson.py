#!/usr/bin/env python3
"""Generate Fig. 7: 18S tree and Coulson-style functional recovery matrices."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo

from plot_style import COLORS, apply_style, figure_mm, normalize_state, panel_label, save_bundle
from tree_plot import draw_tree
from matplotlib.patches import Patch, Wedge


MITO_COLORS = {
    "CI": "#31688E", "CII": "#21918C", "CIII": "#5EAC46", "CIV": "#D9A441",
    "CV": "#E67E22", "TCA": "#756BB1", "PDH": "#C65D7B", "ISC": "#8C6D31",
}
INVASION_COLORS = {
    "Adhesion": "#31688E", "Proteases": "#C65D7B", "Lectin": "#D9A441",
    "Surface": "#5EAC46", "Secreted": "#756BB1",
}


def canonical(value: object) -> str:
    return str(value).strip().replace(" ", "_")


def matrix_lookup(path: Path) -> dict[str, dict[str, object]]:
    data = pd.read_csv(path, sep="\t")
    taxon_col = "taxon" if "taxon" in data.columns else "species"
    return {canonical(row[taxon_col]): row.to_dict() for _, row in data.iterrows()}


def draw_sector_pie(ax, x: float, y: float, values: list[object], recovered_color: str, radius: float = 0.34) -> None:
    width = 360 / len(values)
    for index, value in enumerate(values):
        state = normalize_state(value)
        color = recovered_color if state == "recovered" else COLORS["grey"] if state == "partial" else COLORS["white"]
        start = 90 - (index + 1) * width
        ax.add_patch(Wedge((x, y), radius, start, start + width, facecolor=color, edgecolor="#333333", linewidth=0.35))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree_newick", type=Path)
    parser.add_argument("mito_state_matrix", type=Path)
    parser.add_argument("mito_targets", type=Path)
    parser.add_argument("invasion_state_matrix", type=Path)
    parser.add_argument("invasion_targets", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--taxon-map", type=Path, help="Optional TSV with tree_taxon and matrix_taxon columns")
    parser.add_argument("--allow-missing-taxa", action="store_true")
    parser.add_argument("--highlight", nargs="*", default=["Myxobolus_episquamalis_QD", "Myxobolus_episquamalis_XM", "Myxobolus_pronini_TZ", "Myxobolus_pronini_XY"])
    parser.add_argument("--recovered-label", default="Recovered")
    parser.add_argument("--width-mm", type=float, default=183)
    parser.add_argument("--height-mm", type=float, default=150)
    args = parser.parse_args()

    tree = Phylo.read(args.tree_newick, "newick")
    if args.taxon_map:
        taxon_map = pd.read_csv(args.taxon_map, sep="\t")
        mapping = dict(zip(taxon_map["tree_taxon"].map(canonical), taxon_map["matrix_taxon"].map(canonical)))
        for leaf in tree.get_terminals():
            leaf.name = mapping.get(canonical(leaf.name), canonical(leaf.name))
    mito = matrix_lookup(args.mito_state_matrix)
    invasion = matrix_lookup(args.invasion_state_matrix)
    mito_targets = pd.read_csv(args.mito_targets, sep="\t").sort_values("overall_order")
    invasion_targets = pd.read_csv(args.invasion_targets, sep="\t")
    modules = list(dict.fromkeys(mito_targets["module"]))
    categories = list(dict.fromkeys(invasion_targets["module"]))
    tree_taxa = {canonical(leaf.name) for leaf in tree.get_terminals()}
    missing_mito = sorted(tree_taxa - set(mito))
    missing_invasion = sorted(tree_taxa - set(invasion))
    if (missing_mito or missing_invasion) and not args.allow_missing_taxa:
        details = []
        if missing_mito:
            details.append("mitochondrial matrix: " + ", ".join(missing_mito))
        if missing_invasion:
            details.append("invasion matrix: " + ", ".join(missing_invasion))
        parser.error("Tree taxa missing from " + "; ".join(details) + ". Supply --taxon-map or use --allow-missing-taxa only for diagnostics.")

    apply_style(6.3)
    fig = figure_mm(args.width_mm, args.height_mm)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.8, 2.3, 1.45], wspace=0.08)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_mito = fig.add_subplot(gs[0, 1])
    ax_inv = fig.add_subplot(gs[0, 2])
    leaf_y = draw_tree(ax_tree, tree, stars={canonical(item) for item in args.highlight}, label_size=5.5, show_support=True)
    panel_label(ax_tree, "a", x=-0.04)

    for x, module in enumerate(modules):
        genes = mito_targets.loc[mito_targets["module"] == module, "gene"].tolist()
        for taxon, y in leaf_y.items():
            row = mito.get(canonical(taxon))
            values = [row.get(gene, "Not recovered") if row else "Not recovered" for gene in genes]
            draw_sector_pie(ax_mito, x, y, values, MITO_COLORS.get(module, COLORS["blue"]))
    ax_mito.set_xlim(-0.55, len(modules) - 0.45)
    ax_mito.set_xticks(range(len(modules)), modules, rotation=45, ha="right", fontweight="bold")
    ax_mito.xaxis.tick_top()
    ax_mito.set_ylim(ax_tree.get_ylim())
    ax_mito.set_yticks([])
    ax_mito.spines[:].set_visible(False)
    ax_mito.tick_params(length=0)
    ax_mito.set_title("Nuclear-encoded mitochondrial functions", fontsize=7, pad=24)
    panel_label(ax_mito, "b", x=-0.03)

    for x, category in enumerate(categories):
        sectors = invasion_targets.loc[invasion_targets["module"] == category, "sector"].tolist()
        columns = [f"{category}:{sector}" for sector in sectors]
        for taxon, y in leaf_y.items():
            row = invasion.get(canonical(taxon))
            values = [row.get(column, "Not recovered") if row else "Not recovered" for column in columns]
            draw_sector_pie(ax_inv, x, y, values, INVASION_COLORS.get(category, COLORS["purple"]))
    ax_inv.set_xlim(-0.55, len(categories) - 0.45)
    ax_inv.set_xticks(range(len(categories)), categories, rotation=45, ha="right", fontweight="bold")
    ax_inv.xaxis.tick_top()
    ax_inv.set_ylim(ax_tree.get_ylim())
    ax_inv.set_yticks([])
    ax_inv.spines[:].set_visible(False)
    ax_inv.tick_params(length=0)
    ax_inv.set_title("Invasion and host-interaction categories", fontsize=7, pad=24)
    panel_label(ax_inv, "c", x=-0.03)

    state_handles = [
        Patch(facecolor=COLORS["dark_grey"], edgecolor="#333333", label=args.recovered_label),
        Patch(facecolor=COLORS["grey"], edgecolor="#333333", label="Partially retained"),
        Patch(facecolor=COLORS["white"], edgecolor="#333333", label="Not recovered"),
    ]
    fig.legend(handles=state_handles, loc="lower center", ncol=3, bbox_to_anchor=(0.50, 0.015), fontsize=6)
    fig.text(0.50, 0.002, "Not recovered indicates no qualifying evidence under the stated search criteria; it does not prove biological absence.", ha="center", fontsize=5.5)
    fig.subplots_adjust(left=0.03, right=0.995, top=0.86, bottom=0.08)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
