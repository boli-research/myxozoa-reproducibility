#!/usr/bin/env python3
"""Rebuild Supplementary Fig. S1 from the original PowerPoint media.

The script preserves the original photographs and microscopy fields. It only
performs panel-specific cropping, two documented 90-degree rotations, and
layout/label rendering. No local retouching, contrast normalization, object
removal, or scale-bar reconstruction is performed.
"""

from __future__ import annotations

import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from PIL import Image


ROOT = Path(__file__).resolve().parent
SOURCE_PPTX = ROOT / "source" / "FigS1_original.pptx"
MEDIA_DIR = ROOT / "work" / "media"
OUTPUT_DIR = ROOT / "output"

FIG_WIDTH_MM = 183.0
FIG_HEIGHT_MM = 223.0
MM_PER_INCH = 25.4


@dataclass(frozen=True)
class Panel:
    label: str
    filename: str
    rotation: str | None = None
    focal_x: float = 0.5
    focal_y: float = 0.5


ROWS = [
    {
        "species": "Myxobolus episquamalis",
        "isolate": "QD",
        "panels": [
            Panel("a", "image1.jpeg", focal_y=0.50),
            Panel("b", "image2.tiff"),
        ],
    },
    {
        "species": "Myxobolus episquamalis",
        "isolate": "XM",
        "panels": [
            Panel("c", "image4.jpeg", focal_x=0.50, focal_y=0.51),
            Panel("d", "image5.jpeg", focal_x=0.42, focal_y=0.50),
            Panel("e", "image3.tiff"),
        ],
    },
    {
        "species": "Myxobolus pronini",
        "isolate": "XY",
        "panels": [
            Panel("f", "image6.jpeg"),
            Panel("g", "image7.jpeg", focal_x=0.51, focal_y=0.50),
            Panel("h", "image8.tiff"),
        ],
    },
    {
        "species": "Myxobolus pronini",
        "isolate": "TZ",
        "panels": [
            Panel("i", "image9.jpeg", rotation="cw"),
            Panel("j", "image10.jpeg", rotation="ccw"),
            Panel("k", "image11.tiff"),
        ],
    },
]


def configure_matplotlib() -> None:
    """Set conservative SCI-compatible typography and export settings."""
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.45,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def extract_media() -> None:
    """Extract the untouched embedded media from the source PPTX."""
    if not SOURCE_PPTX.exists():
        raise FileNotFoundError(f"Missing source PPTX: {SOURCE_PPTX}")
    MEDIA_DIR.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(SOURCE_PPTX) as archive:
        for member in archive.namelist():
            if not member.startswith("ppt/media/") or member.endswith("/"):
                continue
            target = MEDIA_DIR / Path(member).name
            with archive.open(member) as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)


def load_panel_image(panel: Panel) -> Image.Image:
    """Load one source image and apply only the documented orientation fix."""
    path = MEDIA_DIR / panel.filename
    image = Image.open(path).convert("RGB")
    if panel.rotation == "cw":
        image = image.transpose(Image.Transpose.ROTATE_270)
    elif panel.rotation == "ccw":
        image = image.transpose(Image.Transpose.ROTATE_90)
    return image


def crop_to_ratio(
    image: Image.Image,
    target_ratio: float,
    focal_x: float = 0.5,
    focal_y: float = 0.5,
) -> Image.Image:
    """Crop to a panel ratio around a documented focal point without rescaling."""
    width, height = image.size
    source_ratio = width / height
    if abs(source_ratio - target_ratio) < 1e-6:
        return image

    if source_ratio > target_ratio:
        crop_width = int(round(height * target_ratio))
        crop_height = height
    else:
        crop_width = width
        crop_height = int(round(width / target_ratio))

    center_x = focal_x * width
    center_y = focal_y * height
    left = int(round(center_x - crop_width / 2))
    top = int(round(center_y - crop_height / 2))
    left = max(0, min(left, width - crop_width))
    top = max(0, min(top, height - crop_height))
    return image.crop((left, top, left + crop_width, top + crop_height))


def add_row_header(fig: plt.Figure, x: float, y: float, species: str, isolate: str) -> None:
    """Draw an italic taxon name followed by a non-italic isolate code."""
    italic = TextArea(
        species,
        textprops={"family": "Arial", "style": "italic", "size": 7.2, "color": "#111111"},
    )
    regular = TextArea(
        f" {isolate}",
        textprops={"family": "Arial", "style": "normal", "weight": "bold", "size": 7.2, "color": "#111111"},
    )
    packed = HPacker(children=[italic, regular], align="center", pad=0, sep=0)
    fig.add_artist(
        AnnotationBbox(
            packed,
            (x, y),
            xycoords="figure fraction",
            box_alignment=(0.0, 0.5),
            frameon=False,
            pad=0,
        )
    )


def add_panel(
    fig: plt.Figure,
    panel: Panel,
    x_mm: float,
    y_mm: float,
    width_mm: float,
    height_mm: float,
) -> None:
    """Place one cropped image with a consistent border and panel label."""
    x = x_mm / FIG_WIDTH_MM
    y = y_mm / FIG_HEIGHT_MM
    width = width_mm / FIG_WIDTH_MM
    height = height_mm / FIG_HEIGHT_MM
    ax = fig.add_axes([x, y, width, height])

    target_ratio = width_mm / height_mm
    image = crop_to_ratio(
        load_panel_image(panel),
        target_ratio=target_ratio,
        focal_x=panel.focal_x,
        focal_y=panel.focal_y,
    )
    ax.imshow(image, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.5, image.width - 0.5)
    ax.set_ylim(image.height - 0.5, -0.5)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.45)
        spine.set_color("#202020")

    ax.text(
        0.014,
        0.982,
        panel.label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        family="Arial",
        fontsize=8.2,
        fontweight="bold",
        color="black",
        bbox={
            "boxstyle": "square,pad=0.13",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.88,
        },
        zorder=5,
    )


def build_figure() -> plt.Figure:
    """Assemble the four isolate groups into one full-page supplementary plate."""
    fig = plt.figure(
        figsize=(FIG_WIDTH_MM / MM_PER_INCH, FIG_HEIGHT_MM / MM_PER_INCH),
        facecolor="white",
    )

    left_margin = 4.0
    right_margin = 4.0
    top_margin = 4.0
    header_height = 4.5
    row_gap = 2.2
    panel_gap = 2.5
    image_heights = [52.2, 46.0, 46.0, 46.0]
    content_width = FIG_WIDTH_MM - left_margin - right_margin

    y_top = FIG_HEIGHT_MM - top_margin
    for row_index, row in enumerate(ROWS):
        header_center_y = y_top - header_height / 2
        add_row_header(
            fig,
            left_margin / FIG_WIDTH_MM,
            header_center_y / FIG_HEIGHT_MM,
            row["species"],
            row["isolate"],
        )

        image_height = image_heights[row_index]
        image_y = y_top - header_height - image_height
        panels = row["panels"]

        if len(panels) == 2:
            # Keep the microscopy panel near its native 4:3 ratio; allocate the
            # remaining width to the whole-host view.
            second_width = 70.0
            first_width = content_width - panel_gap - second_width
            widths = [first_width, second_width]
        elif row_index == 1:
            # The XM host image is intrinsically wide. Preserve the complete
            # fish in panel c, keep the cyst close-up compact, and leave the
            # microscopy field close to its native 4:3 aspect ratio.
            widths = [68.0, 42.0, content_width - panel_gap * 2 - 110.0]
        else:
            width = (content_width - panel_gap * 2) / 3
            widths = [width, width, width]

        x = left_margin
        for panel, width in zip(panels, widths):
            add_panel(fig, panel, x, image_y, width, image_height)
            x += width + panel_gap

        y_top = image_y - row_gap

    return fig


def export(fig: plt.Figure) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = OUTPUT_DIR / "Supplementary_Figure_S1"

    fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches=None)
    fig.savefig(stem.with_suffix(".svg"), format="svg", bbox_inches=None)
    fig.savefig(stem.with_suffix(".png"), format="png", dpi=300, bbox_inches=None)
    fig.savefig(
        stem.with_suffix(".tiff"),
        format="tiff",
        dpi=600,
        bbox_inches=None,
        pil_kwargs={"compression": "tiff_lzw"},
    )


def main() -> None:
    configure_matplotlib()
    extract_media()
    figure = build_figure()
    export(figure)
    plt.close(figure)
    print(f"Exports written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
