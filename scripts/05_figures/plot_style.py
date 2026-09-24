from __future__ import annotations

import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "myxobolus-matplotlib-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MM_TO_INCH = 1 / 25.4
COLORS = {
    "blue": "#31688E",
    "teal": "#21918C",
    "green": "#5EAC46",
    "gold": "#D9A441",
    "orange": "#E67E22",
    "rose": "#C65D7B",
    "purple": "#756BB1",
    "grey": "#A6A6A6",
    "light_grey": "#E6E6E6",
    "dark_grey": "#4D4D4D",
    "black": "#202020",
    "white": "#FFFFFF",
}


def apply_style(font_size: float = 7.0) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": font_size,
        "axes.linewidth": 0.8,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    })


def figure_mm(width_mm: float, height_mm: float, **kwargs):
    return plt.figure(figsize=(width_mm * MM_TO_INCH, height_mm * MM_TO_INCH), **kwargs)


def panel_label(ax, label: str, x: float = -0.08, y: float = 1.02) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontweight="bold", fontsize=8, ha="left", va="bottom")


def save_bundle(fig, output_prefix: str | Path, dpi: int = 600) -> list[Path]:
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix, options in (
        ("svg", {}),
        ("pdf", {}),
        ("tiff", {"dpi": dpi, "pil_kwargs": {"compression": "tiff_lzw"}}),
        ("png", {"dpi": 300}),
    ):
        path = prefix.with_suffix(f".{suffix}")
        fig.savefig(path, bbox_inches="tight", **options)
        outputs.append(path)
    return outputs


def normalize_state(value: object) -> str:
    text = str(value).strip().lower().replace("_", " ")
    if text in {"1", "1.0", "present", "recovered", "complete"}:
        return "recovered"
    if text in {"0.5", "partial", "partially retained", "partially recovered"}:
        return "partial"
    if text in {"0", "0.0", "absent", "missing", "not recovered", "nan", "none", ""}:
        return "missing"
    raise ValueError(f"Unrecognized recovery state: {value!r}")
