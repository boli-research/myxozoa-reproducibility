#!/usr/bin/env python3
"""Generate Fig. 3: mitochondrial tree plus standardized gene-order tracks."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo

from plot_style import COLORS, apply_style, figure_mm, panel_label, save_bundle
from tree_plot import draw_tree
from matplotlib.patches import FancyArrow, Patch


FEATURE_COLORS = {
    "cox1": "#31688E", "cox2": "#4C9ED9", "nad1": "#5EAC46", "nad5": "#91CF60",
    "cytb": "#E67E22", "rrnS": "#756BB1", "rrnL": "#A58CC7", "ORF": "#A6A6A6",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree_newick", type=Path)
    parser.add_argument("features_tsv", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--taxon-map", type=Path, help="Optional TSV with tree_taxon and feature_taxon columns")
    args = parser.parse_args()

    tree = Phylo.read(args.tree_newick, "newick")
    data = pd.read_csv(args.features_tsv, sep="\t")
    if args.taxon_map:
        taxon_map = pd.read_csv(args.taxon_map, sep="\t")
        mapping = dict(zip(taxon_map["tree_taxon"].astype(str), taxon_map["feature_taxon"].astype(str)))
        for leaf in tree.get_terminals():
            leaf.name = mapping.get(str(leaf.name), str(leaf.name))
    tree_taxa = {str(leaf.name).replace(" ", "_") for leaf in tree.get_terminals()}
    feature_taxa = set(data["taxon"].astype(str).str.replace(" ", "_"))
    missing = sorted(tree_taxa - feature_taxa)
    if missing:
        raise ValueError("No mitochondrial features for tree taxa: " + ", ".join(missing))
    apply_style(6.2)
    height = max(90, 7.2 * len(tree.get_terminals()))
    fig = figure_mm(183, height)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 2.0], wspace=0.03)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_order = fig.add_subplot(gs[0, 1])
    leaf_y = draw_tree(ax_tree, tree, label_size=5.5, show_support=True)
    panel_label(ax_tree, "a", x=-0.03)

    max_end = float(data["end"].max())
    for taxon, y in leaf_y.items():
        subset = data[data["taxon"].astype(str).str.replace(" ", "_") == str(taxon).replace(" ", "_")]
        ax_order.plot([0, max_end], [y, y], color="#D0D0D0", linewidth=0.45, zorder=0)
        for _, feature in subset.iterrows():
            start = float(feature["start"])
            length = float(feature["end"] - feature["start"] + 1)
            strand = int(feature.get("strand", 1))
            if strand < 0:
                start += length
                length = -length
            head = min(abs(length) * 0.25, 600)
            patch = FancyArrow(start, y - 0.18, length, 0, width=0.36, head_width=0.36, head_length=head, length_includes_head=True, color=FEATURE_COLORS.get(feature["feature"], COLORS["grey"]), linewidth=0)
            ax_order.add_patch(patch)
    ax_order.set_xlim(0, max_end * 1.02)
    ax_order.set_ylim(ax_tree.get_ylim())
    ax_order.set_yticks([])
    ax_order.set_xlabel("Linearized mitochondrial position (bp)")
    ax_order.spines[["left", "right", "top"]].set_visible(False)
    panel_label(ax_order, "b", x=-0.03)
    handles = [Patch(facecolor=color, label=name) for name, color in FEATURE_COLORS.items()]
    ax_order.legend(handles=handles, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=5.5)
    fig.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.12)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
