#!/usr/bin/env python3
"""Generate Fig. 2 circular mitochondrial genome maps from GenBank files."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np
from Bio import SeqIO

from plot_style import COLORS, apply_style, figure_mm, panel_label, save_bundle


GENE_COLORS = {
    "cox1": "#31688E", "cox2": "#4C9ED9", "nad1": "#5EAC46", "nad5": "#91CF60",
    "cytb": "#E67E22", "rrnS": "#756BB1", "rrnL": "#A58CC7", "ORF": "#A6A6A6",
}


def feature_name(feature) -> str:
    values = []
    for key in ("gene", "product", "label", "locus_tag"):
        values.extend(feature.qualifiers.get(key, []))
    text = " ".join(values).lower()
    aliases = [("cox1", "cox1"), ("cox2", "cox2"), ("nad1", "nad1"), ("nad5", "nad5"), ("cob", "cytb"), ("cytb", "cytb")]
    for token, name in aliases:
        if re.search(rf"\b{token}\b", text):
            return name
    if feature.type == "rRNA":
        return "rrnS" if "small" in text or "12s" in text else "rrnL"
    return "ORF"


def circular_profile(sequence: str, window: int = 1000, step: int = 100):
    sequence = sequence.upper()
    positions, gc, skew = [], [], []
    for start in range(0, len(sequence), step):
        fragment = "".join(sequence[(start + offset) % len(sequence)] for offset in range(window))
        g, c = fragment.count("G"), fragment.count("C")
        valid = sum(fragment.count(base) for base in "ACGT")
        positions.append(start)
        gc.append((g + c) / valid if valid else 0)
        skew.append((g - c) / (g + c) if g + c else 0)
    return np.asarray(positions), np.asarray(gc), np.asarray(skew)


def selected_features(record, min_orf_bp: int = 150, max_overlap_fraction: float = 0.10):
    core = []
    orfs = []
    for feature in record.features:
        if feature.type not in {"CDS", "rRNA"}:
            continue
        name = feature_name(feature)
        if name == "ORF":
            orfs.append((feature, name))
        else:
            core.append((feature, name))
    core_positions = [(feature, set(int(position) for position in feature.location)) for feature, _ in core]
    retained = list(core)
    for feature, name in orfs:
        positions = set(int(position) for position in feature.location)
        if len(positions) <= min_orf_bp:
            continue
        overlap_too_large = False
        for _, reference in core_positions:
            denominator = min(len(positions), len(reference))
            if denominator and len(positions & reference) / denominator > max_overlap_fraction:
                overlap_too_large = True
                break
        if not overlap_too_large:
            retained.append((feature, name))
    return sorted(retained, key=lambda item: int(item[0].location.start))


def draw_record(ax, record, title: str, source_dir: Path | None = None) -> None:
    length = len(record.seq)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.05)
    ax.set_axis_off()
    orf_index = 0
    feature_rows = []
    for feature, name in selected_features(record):
        if name == "ORF":
            orf_index += 1
            display = f"ORF{orf_index}"
        else:
            display = name
        start = int(feature.location.start)
        end = int(feature.location.end)
        theta_start = 2 * np.pi * start / length
        theta_width = 2 * np.pi * max(end - start, 1) / length
        bottom = 0.78 if (feature.location.strand or 1) > 0 else 0.68
        ax.bar(theta_start + theta_width / 2, 0.085, width=theta_width, bottom=bottom, color=GENE_COLORS[name], edgecolor="white", linewidth=0.35)
        if name != "ORF" or end - start >= 300:
            angle = theta_start + theta_width / 2
            rotation = np.degrees(np.pi / 2 - angle)
            if rotation < -90:
                rotation += 180
            if rotation > 90:
                rotation -= 180
            ax.text(angle, 0.91, display, fontsize=4.5, rotation=rotation, rotation_mode="anchor", ha="center", va="center")
        feature_rows.append([title, display, name, start + 1, end, feature.location.strand or 0])

    positions, gc, skew = circular_profile(str(record.seq))
    theta = 2 * np.pi * positions / length
    gc_scaled = 0.55 + (gc - gc.mean()) * 0.9
    ax.plot(theta, gc_scaled, color=COLORS["teal"], linewidth=0.55)
    ax.axhline(0.55, color="#CFCFCF", linewidth=0.35)
    positive = np.clip(skew, 0, None)
    negative = np.clip(skew, None, 0)
    ax.fill_between(theta, 0.38, 0.38 + positive * 0.18, color="#5EAC46", linewidth=0)
    ax.fill_between(theta, 0.38, 0.38 + negative * 0.18, color="#C65D7B", linewidth=0)
    ax.text(0, 0.13, f"{length:,} bp\nGC {100 * (str(record.seq).upper().count('G') + str(record.seq).upper().count('C')) / length:.1f}%", ha="center", va="center", fontsize=6)
    ax.set_title(title.replace("_", " "), fontsize=7, fontstyle="italic", pad=8)

    if source_dir:
        source_dir.mkdir(parents=True, exist_ok=True)
        with (source_dir / f"{title}.features.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["taxon", "display_name", "feature", "start", "end", "strand"])
            writer.writerows(feature_rows)
        with (source_dir / f"{title}.gc_profile.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["position_bp", "gc_content", "gc_skew"])
            writer.writerows(zip(positions + 1, gc, skew))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV with taxon and genbank columns")
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--source-data-dir", type=Path)
    args = parser.parse_args()
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 4:
        parser.error("The final Fig. 2 layout expects exactly four isolates")

    apply_style(6.5)
    fig = figure_mm(183, 100)
    axes = [fig.add_subplot(1, 4, index + 1, projection="polar") for index in range(4)]
    for index, (ax, row) in enumerate(zip(axes, rows)):
        records = list(SeqIO.parse(row["genbank"], "genbank"))
        if len(records) != 1:
            raise ValueError(f"Expected one GenBank record: {row['genbank']}")
        draw_record(ax, records[0], row["taxon"], args.source_data_dir)
        panel_label(ax, chr(ord("a") + index), x=-0.05, y=1.02)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.04, wspace=0.15)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
