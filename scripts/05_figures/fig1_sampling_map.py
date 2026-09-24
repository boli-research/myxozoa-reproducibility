#!/usr/bin/env python3
"""Generate Fig. 1 sampling maps from public coordinates and GeoJSON boundaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from PIL import Image

from plot_style import apply_style, figure_mm, panel_label, save_bundle
from matplotlib.offsetbox import AnnotationBbox, OffsetImage


def polygons(geometry: dict):
    if geometry["type"] == "Polygon":
        yield geometry["coordinates"]
    elif geometry["type"] == "MultiPolygon":
        yield from geometry["coordinates"]


def load_polygons(path: Path) -> list[list]:
    data = json.loads(path.read_text(encoding="utf-8"))
    features = data.get("features", [data])
    output = []
    for feature in features:
        output.extend(polygons(feature.get("geometry", feature)))
    return output


def draw_boundaries(ax, shapes: list[list], color: str = "#666666", linewidth: float = 0.45) -> None:
    for polygon in shapes:
        exterior = polygon[0]
        x, y = zip(*[(point[0], point[1]) for point in exterior])
        ax.plot(x, y, color=color, linewidth=linewidth, zorder=1)


def bounds(shapes: list[list]) -> tuple[float, float, float, float]:
    points = [point for polygon in shapes for ring in polygon for point in ring]
    x = [point[0] for point in points]
    y = [point[1] for point in points]
    return min(x), max(x), min(y), max(y)


def add_optional_image(ax, path: str, xy: tuple[float, float], zoom: float = 0.13) -> None:
    if not path or str(path).lower() in {"nan", "na"} or not Path(path).is_file():
        return
    image = Image.open(path).convert("RGBA")
    ax.add_artist(AnnotationBbox(OffsetImage(image, zoom=zoom), xy, frameon=False, box_alignment=(0.5, 0.5)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sites_tsv", type=Path)
    parser.add_argument("china_geojson", type=Path)
    parser.add_argument("province_dir", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--artwork-manifest", type=Path, help="Optional TSV: isolate, host_image, spore_image")
    args = parser.parse_args()

    sites = pd.read_csv(args.sites_tsv, sep="\t")
    artwork = {}
    if args.artwork_manifest:
        art = pd.read_csv(args.artwork_manifest, sep="\t")
        artwork = {row["isolate"]: row.to_dict() for _, row in art.iterrows()}

    apply_style(6.5)
    fig = figure_mm(183, 112)
    gs = fig.add_gridspec(2, 4, height_ratios=[1.05, 1], hspace=0.16, wspace=0.12)
    ax_china = fig.add_subplot(gs[0, :])
    local_axes = [fig.add_subplot(gs[1, index]) for index in range(4)]

    china = load_polygons(args.china_geojson)
    draw_boundaries(ax_china, china, linewidth=0.35)
    ax_china.scatter(sites["longitude_E"], sites["latitude_N"], marker="v", s=22, color="#C43C39", edgecolor="white", linewidth=0.4, zorder=3)
    for _, site in sites.iterrows():
        ax_china.text(site["longitude_E"] + 0.7, site["latitude_N"] + 0.5, site["isolate"], fontsize=6, fontweight="bold")
    ax_china.set_xlim(73, 135)
    ax_china.set_ylim(18, 54)
    ax_china.set_aspect("equal")
    ax_china.axis("off")
    panel_label(ax_china, "a", x=0.005, y=0.98)

    for ax, (_, site) in zip(local_axes, sites.iterrows()):
        province_file = args.province_dir / f"{str(site['province']).lower()}.geojson"
        if not province_file.is_file():
            raise FileNotFoundError(f"Province boundary not found: {province_file}")
        shapes = load_polygons(province_file)
        draw_boundaries(ax, shapes, linewidth=0.55)
        x0, x1, y0, y1 = bounds(shapes)
        xpad = max((x1 - x0) * 0.08, 0.2)
        ypad = max((y1 - y0) * 0.08, 0.2)
        ax.set_xlim(x0 - xpad, x1 + xpad)
        ax.set_ylim(y0 - ypad, y1 + ypad)
        ax.scatter(site["longitude_E"], site["latitude_N"], marker="v", s=30, color="#C43C39", edgecolor="white", linewidth=0.5, zorder=3)
        ax.text(0.03, 0.97, f"{site['city']}, {site['province']}\n{site['isolate']} · {site['habitat']}", transform=ax.transAxes, va="top", fontsize=5.5)
        ax.text(0.03, 0.04, f"{site['host']} · {site['parasitic_site']}", transform=ax.transAxes, fontsize=5.2, fontstyle="italic")
        record = artwork.get(site["isolate"], {})
        add_optional_image(ax, record.get("host_image", ""), (x0 + 0.25 * (x1 - x0), y0 + 0.25 * (y1 - y0)))
        add_optional_image(ax, record.get("spore_image", ""), (x0 + 0.75 * (x1 - x0), y0 + 0.25 * (y1 - y0)))
        ax.set_aspect("equal")
        ax.axis("off")
        panel_label(ax, str(site["panel"]).lower(), x=0.01, y=0.98)

    fig.subplots_adjust(left=0.02, right=0.99, top=0.98, bottom=0.02)
    save_bundle(fig, args.output_prefix)


if __name__ == "__main__":
    main()
