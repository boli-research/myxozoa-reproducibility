#!/usr/bin/env python3
"""Generate Fig. 5: assembly statistics and BUSCO completeness."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from plot_style import COLORS, apply_style, figure_mm, panel_label, save_bundle
from matplotlib.patches import Patch


GROUP_COLORS = {
    "Myxobolidae": "#31688E",
    "Kudoidae": "#7A6FAC",
    "Ceratomyxidae": "#5EAC46",
    "Other Myxozoa": "#A6A6A6",
    "new": "#C65D7B",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_tsv", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--width-mm", type=float, default=183)
    parser.add_argument("--height-mm", type=float, default=105)
    args = parser.parse_args()

    data = pd.read_csv(args.source_tsv, sep="\t")
    required = {
        "taxon", "group", "assembly_size_mb", "contigs", "n50_bp", "gc_percent",
        "busco_single", "busco_duplicated", "busco_fragmented", "busco_missing",
    }
    missing = required - set(data.columns)
    if missing:
        parser.error("Missing columns: " + ", ".join(sorted(missing)))
    data = data.sort_values(["group", "taxon"], kind="stable").reset_index(drop=True)
    y = np.arange(len(data))[::-1]
    colors = [GROUP_COLORS.get(group, COLORS["grey"]) for group in data["group"]]

    apply_style(6.5)
    fig = figure_mm(args.width_mm, args.height_mm)
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.1], hspace=0.36, wspace=0.42)
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    ax_busco = fig.add_subplot(gs[1, :])

    metrics = [
        ("assembly_size_mb", "Assembly size (Mb)", False),
        ("n50_bp", "Contig N50 (bp)", True),
        ("contigs", "Number of contigs", True),
        ("gc_percent", "GC content (%)", False),
    ]
    for index, (ax, (column, xlabel, log_scale)) in enumerate(zip(axes, metrics)):
        ax.barh(y, data[column], color=colors, edgecolor="white", height=0.72)
        ax.set_xlabel(xlabel)
        if log_scale:
            ax.set_xscale("log")
        if index == 0:
            ax.set_yticks(y, [name.replace("_", " ") for name in data["taxon"]], fontsize=5.5, fontstyle="italic")
        else:
            ax.set_yticks([])
        ax.grid(axis="x", color="#E8E8E8", linewidth=0.5, zorder=0)
        panel_label(ax, chr(ord("a") + index))

    components = [
        ("busco_single", "Complete, single-copy", "#3B7EA1"),
        ("busco_duplicated", "Complete, duplicated", "#7DB7D5"),
        ("busco_fragmented", "Fragmented", "#D9A441"),
        ("busco_missing", "Missing", "#D9D9D9"),
    ]
    left = np.zeros(len(data))
    for column, label, color in components:
        values = data[column].to_numpy(float)
        ax_busco.barh(y, values, left=left, color=color, edgecolor="white", height=0.72, label=label)
        left += values
    ax_busco.set_yticks(y, [name.replace("_", " ") for name in data["taxon"]], fontsize=5.5, fontstyle="italic")
    ax_busco.set_xlim(0, 100)
    ax_busco.set_xlabel("BUSCO genes (%)")
    ax_busco.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6)
    panel_label(ax_busco, "e", x=-0.03)

    group_handles = [Patch(facecolor=color, label=group) for group, color in GROUP_COLORS.items() if group in set(data["group"])]
    if group_handles:
        axes[-1].legend(handles=group_handles, loc="best", fontsize=5.5)
    fig.subplots_adjust(left=0.19, right=0.99, top=0.96, bottom=0.18)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
