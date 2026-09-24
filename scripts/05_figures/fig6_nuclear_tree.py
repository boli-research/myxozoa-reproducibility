#!/usr/bin/env python3
"""Generate Fig. 6 nuclear phylogeny with lineage colors and new-isolate markers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo

from plot_style import apply_style, figure_mm, panel_label, save_bundle
from tree_plot import draw_tree
from matplotlib.patches import Patch


DEFAULT_COLORS = ["#31688E", "#5EAC46", "#D9A441", "#C65D7B", "#756BB1", "#4C9ED9", "#777777"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree_newick", type=Path)
    parser.add_argument("taxon_metadata", type=Path, help="TSV: taxon, lineage, new_isolate")
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()

    tree = Phylo.read(args.tree_newick, "newick")
    metadata = pd.read_csv(args.taxon_metadata, sep="\t").fillna("")
    lineages = list(dict.fromkeys(metadata["lineage"]))
    lineage_colors = {lineage: DEFAULT_COLORS[index % len(DEFAULT_COLORS)] for index, lineage in enumerate(lineages)}
    metadata["key"] = metadata["taxon"].astype(str).str.replace(" ", "_")
    tree_taxa = {str(leaf.name).replace(" ", "_") for leaf in tree.get_terminals()}
    missing = sorted(tree_taxa - set(metadata["key"]))
    if missing:
        raise ValueError("No lineage metadata for tree taxa: " + ", ".join(missing))
    label_colors = dict(zip(metadata["key"], metadata["lineage"].map(lineage_colors)))
    stars = set(metadata.loc[metadata["new_isolate"].astype(str).str.lower().isin(["1", "true", "yes"]), "key"])

    apply_style(6.0)
    height = max(130, 3.6 * len(tree.get_terminals()))
    fig = figure_mm(125, height)
    ax = fig.add_subplot(111)
    draw_tree(ax, tree, label_colors=label_colors, stars=stars, label_size=5.3, show_support=True)
    panel_label(ax, "a", x=-0.02)
    handles = [Patch(facecolor=color, label=lineage) for lineage, color in lineage_colors.items()]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, -0.08), ncol=min(4, len(handles)), fontsize=5.5)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.98, bottom=0.10)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
