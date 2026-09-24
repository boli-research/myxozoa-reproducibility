#!/usr/bin/env python3
"""Generate Fig. 4: cnidarian mitochondrial phylogeny and chromosome architecture."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo

from plot_style import COLORS, apply_style, figure_mm, panel_label, save_bundle
from tree_plot import draw_tree
from matplotlib.patches import Circle, Patch


LINEAGE_COLORS = {
    "Hexacorallia": "#31688E", "Octocorallia": "#4C9ED9", "Staurozoa": "#756BB1",
    "Scyphozoa": "#C65D7B", "Hydrozoa": "#5EAC46", "Polypodiozoa": "#D9A441",
    "Myxozoa": "#E67E22", "Porifera": "#777777",
}


def architecture_marker(ax, x: float, y: float, architecture: str) -> None:
    architecture = architecture.lower().replace(" ", "_")
    if architecture in {"single_circular", "circular"}:
        ax.add_patch(Circle((x, y), 0.20, facecolor="white", edgecolor=COLORS["dark_grey"], linewidth=0.8))
    elif architecture in {"multiple_circular", "multipartite_circular"}:
        for offset in (-0.22, 0, 0.22):
            ax.add_patch(Circle((x + offset, y), 0.11, facecolor="white", edgecolor=COLORS["dark_grey"], linewidth=0.7))
    elif architecture in {"single_linear", "linear"}:
        ax.plot([x - 0.30, x + 0.30], [y, y], color=COLORS["dark_grey"], linewidth=1.3)
    elif architecture in {"multiple_linear", "multipartite_linear"}:
        for offset in (-0.12, 0.12):
            ax.plot([x - 0.30, x + 0.30], [y + offset, y + offset], color=COLORS["dark_grey"], linewidth=1.0)
    elif architecture in {"absent", "mitochondrial_genome_absent"}:
        ax.text(x, y, "×", color="#C43C39", ha="center", va="center", fontsize=8, fontweight="bold")
    else:
        ax.text(x, y, "?", color=COLORS["grey"], ha="center", va="center", fontsize=7)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree_newick", type=Path)
    parser.add_argument("taxon_metadata", type=Path, help="TSV: taxon, lineage, architecture, new_isolate")
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()

    tree = Phylo.read(args.tree_newick, "newick")
    metadata = pd.read_csv(args.taxon_metadata, sep="\t").fillna("")
    metadata["key"] = metadata["taxon"].astype(str).str.replace(" ", "_")
    lookup = metadata.set_index("key").to_dict("index")
    tree_taxa = {str(leaf.name).replace(" ", "_") for leaf in tree.get_terminals()}
    missing = sorted(tree_taxa - set(lookup))
    if missing:
        raise ValueError("No architecture metadata for tree taxa: " + ", ".join(missing))
    label_colors = {key: LINEAGE_COLORS.get(row["lineage"], "black") for key, row in lookup.items()}
    stars = {key for key, row in lookup.items() if str(row.get("new_isolate", "")).lower() in {"1", "true", "yes"}}

    apply_style(6.1)
    height = max(130, 4.0 * len(tree.get_terminals()))
    fig = figure_mm(183, height)
    gs = fig.add_gridspec(1, 2, width_ratios=[3.5, 0.7], wspace=0.02)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax_arch = fig.add_subplot(gs[0, 1])
    leaf_y = draw_tree(ax_tree, tree, label_colors=label_colors, stars=stars, label_size=5.2, show_support=True)
    panel_label(ax_tree, "a", x=-0.02)
    for taxon, y in leaf_y.items():
        row = lookup.get(str(taxon).replace(" ", "_"), {})
        architecture_marker(ax_arch, 0, y, str(row.get("architecture", "unknown")))
    ax_arch.set_ylim(ax_tree.get_ylim())
    ax_arch.set_xlim(-0.6, 0.6)
    ax_arch.set_xticks([0], ["mtDNA\narchitecture"], fontsize=6)
    ax_arch.xaxis.tick_top()
    ax_arch.set_yticks([])
    ax_arch.spines[:].set_visible(False)
    panel_label(ax_arch, "b", x=-0.12)
    handles = [Patch(facecolor=color, label=lineage) for lineage, color in LINEAGE_COLORS.items() if lineage in set(metadata["lineage"])]
    fig.legend(handles=handles, loc="lower center", ncol=min(5, len(handles)), fontsize=5.5)
    fig.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.08)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
