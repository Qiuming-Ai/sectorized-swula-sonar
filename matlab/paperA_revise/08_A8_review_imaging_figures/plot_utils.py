"""Shared publication plotting style for the Paper-A review figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


METHOD_ORDER = ["dense", "cdm", "fdm", "random", "nested", "coprime"]
METHOD_LABELS = {
    "dense": "Dense ULA",
    "cdm": "CDM SW-ULA",
    "fdm": "FDM SW-ULA",
    "random": "Random sparse",
    "nested": "Nested sparse",
    "coprime": "Coprime sparse",
}

# cividis is perceptually uniform and colour-vision-deficiency safe.
IMAGE_CMAP = "cividis"
TARGET_COLOR = "white"
ZOOM_COLOR = "#D55E00"  # Okabe-Ito vermilion.


def apply_paper_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.8,
            "ytick.major.size": 2.8,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, output_base: Path, dpi: int = 600) -> None:
    """Save the requested high-resolution PNG plus a vector-text PDF copy."""
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=dpi, facecolor="white")
    fig.savefig(output_base.with_suffix(".pdf"), dpi=300, facecolor="white")

