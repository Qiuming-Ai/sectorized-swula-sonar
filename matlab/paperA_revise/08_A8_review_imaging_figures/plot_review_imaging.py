"""Render reviewer-ready Paper-A angle-range images from MATLAB exports.

Outputs
-------
* Fig16_imaging_comparison_full_and_zoom.png/pdf
* Fig16_imaging_comparison_full.png/pdf
* Fig16_imaging_comparison_zoom.png/pdf
* Fig17_low_snr_comparison.png/pdf
* individual/<method>_snr_m15_full_and_zoom.png/pdf

Every image uses an angle-range coordinate system, a shared dynamic range,
physical units, a labelled colour bar, ground-truth target markers, and an
explicitly dimensioned zoom window.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.io import loadmat

from plot_utils import (
    IMAGE_CMAP,
    METHOD_LABELS,
    METHOD_ORDER,
    TARGET_COLOR,
    ZOOM_COLOR,
    apply_paper_style,
    save_figure,
)


FULL_BOUNDS = (-50.0, 50.0, 4.80, 5.70)
ZOOM_BOUNDS = (-8.0, 8.0, 4.96, 5.32)
DEFAULT_DYNAMIC_RANGE_DB = 40.0


@dataclass(frozen=True)
class ImageData:
    suite: str
    method_key: str
    method_name: str
    snr_db: float
    angle_deg: np.ndarray
    range_m: np.ndarray
    image_db: np.ndarray
    target_angle_deg: np.ndarray
    target_range_m: np.ndarray
    source_script: str


def _mat_string(value: np.ndarray | str) -> str:
    if isinstance(value, str):
        return value
    array = np.asarray(value).squeeze()
    if array.ndim == 0:
        return str(array.item())
    if array.dtype.kind in {"U", "S"}:
        return "".join(array.tolist()).strip()
    return str(array.tolist())


def load_image(path: Path, dynamic_range_db: float) -> ImageData:
    raw = loadmat(path, squeeze_me=True)
    angle = np.asarray(raw["angle_deg"], dtype=float).reshape(-1)
    ranges = np.asarray(raw["range_m"], dtype=float).reshape(-1)
    image = np.abs(np.asarray(raw["image_linear"]))
    if image.shape == (angle.size, ranges.size):
        image = image.T
    if image.shape != (ranges.size, angle.size):
        raise ValueError(
            f"{path.name}: image {image.shape}, expected {(ranges.size, angle.size)}"
        )

    finite_max = np.nanmax(image)
    if not np.isfinite(finite_max) or finite_max <= 0:
        raise ValueError(f"{path.name}: image has no positive finite samples")
    image_db = 20.0 * np.log10(np.maximum(image / finite_max, np.finfo(float).tiny))
    image_db = np.clip(image_db, -dynamic_range_db, 0.0)

    return ImageData(
        suite=_mat_string(raw["suite"]),
        method_key=_mat_string(raw["method_key"]),
        method_name=_mat_string(raw["method_name"]),
        snr_db=float(np.asarray(raw["snr_db"]).squeeze()),
        angle_deg=angle,
        range_m=ranges,
        image_db=image_db,
        target_angle_deg=np.asarray(raw["target_angle_deg"], dtype=float).reshape(-1),
        target_range_m=np.asarray(raw["target_range_m"], dtype=float).reshape(-1),
        source_script=_mat_string(raw["source_script"]),
    )


def load_dataset(data_dir: Path, dynamic_range_db: float) -> dict[tuple[str, str, int], ImageData]:
    dataset: dict[tuple[str, str, int], ImageData] = {}
    for path in sorted(data_dir.glob("*.mat")):
        item = load_image(path, dynamic_range_db)
        key = (item.suite, item.method_key, int(round(item.snr_db)))
        if key in dataset:
            raise ValueError(f"duplicate data key {key}: {path}")
        dataset[key] = item
    if not dataset:
        raise FileNotFoundError(f"No MATLAB exports found in {data_dir}")
    return dataset


def plot_image(
    ax: plt.Axes,
    item: ImageData,
    bounds: tuple[float, float, float, float],
    dynamic_range_db: float,
    *,
    show_targets: bool = True,
    show_zoom_box: bool = False,
) -> object:
    x0, x1, y0, y1 = bounds
    mesh = ax.pcolormesh(
        item.angle_deg,
        item.range_m,
        item.image_db,
        shading="auto",
        cmap=IMAGE_CMAP,
        vmin=-dynamic_range_db,
        vmax=0.0,
        rasterized=True,
    )
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_xticks(np.arange(np.ceil(x0 / 10) * 10, x1 + 0.1, 10))
    y_step = 0.2 if (y1 - y0) > 0.5 else 0.1
    ax.set_yticks(np.arange(np.ceil(y0 / y_step) * y_step, y1 + 1e-9, y_step))
    ax.grid(False)

    if show_targets:
        mask = (
            (item.target_angle_deg >= x0)
            & (item.target_angle_deg <= x1)
            & (item.target_range_m >= y0)
            & (item.target_range_m <= y1)
        )
        ax.scatter(
            item.target_angle_deg[mask],
            item.target_range_m[mask],
            marker="+",
            s=22,
            linewidths=0.85,
            c=TARGET_COLOR,
            label="Ground-truth target",
            zorder=4,
        )

    if show_zoom_box:
        zx0, zx1, zy0, zy1 = ZOOM_BOUNDS
        ax.add_patch(
            Rectangle(
                (zx0, zy0),
                zx1 - zx0,
                zy1 - zy0,
                fill=False,
                edgecolor=ZOOM_COLOR,
                linewidth=1.1,
                linestyle="--",
                zorder=5,
            )
        )
        ax.text(
            zx1 - 0.6,
            zy0 + 0.018,
            "16° × 0.36 m",
            color=ZOOM_COLOR,
            fontsize=6.8,
            ha="right",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.8},
            zorder=6,
        )
    return mesh


def add_shared_colorbar(fig: plt.Figure, mesh: object, axes: np.ndarray, dynamic_range_db: float) -> None:
    colorbar = fig.colorbar(mesh, ax=np.asarray(axes).ravel().tolist(), pad=0.012, fraction=0.025)
    colorbar.set_label("Normalized amplitude (dB)")
    colorbar.set_ticks(np.arange(-dynamic_range_db, 0.1, 10))


def main_composite(dataset: dict, output_dir: Path, dynamic_range_db: float) -> None:
    fig, axes = plt.subplots(6, 2, figsize=(7.2, 10.0), constrained_layout=True)
    panel_letters = iter("abcdefghijkl")
    mesh = None
    for row, method in enumerate(METHOD_ORDER):
        item = dataset[("main", method, -15)]
        left, right = axes[row]
        mesh = plot_image(left, item, FULL_BOUNDS, dynamic_range_db, show_zoom_box=True)
        plot_image(right, item, ZOOM_BOUNDS, dynamic_range_db)
        left.set_title(f"({next(panel_letters)}) {METHOD_LABELS[method]}: full FOV")
        right.set_title(f"({next(panel_letters)}) {METHOD_LABELS[method]}: zoom")
        left.set_ylabel("Range (m)")
        right.set_ylabel("Range (m)")
        if row == 5:
            left.set_xlabel("Angle (deg)")
            right.set_xlabel("Angle (deg)")
    assert mesh is not None
    add_shared_colorbar(fig, mesh, axes, dynamic_range_db)
    fig.text(
        0.5,
        1.002,
        f"Input SNR = −15 dB; FOV = 100°; display dynamic range = {dynamic_range_db:g} dB; + denotes ground truth",
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    save_figure(fig, output_dir / "Fig16_imaging_comparison_full_and_zoom")
    plt.close(fig)


def three_method_composite(dataset: dict, output_dir: Path, dynamic_range_db: float) -> None:
    """Main-paper comparison restricted to Dense/CDM/FDM."""
    methods = ["dense", "cdm", "fdm"]
    fig, axes = plt.subplots(3, 2, figsize=(7.2, 6.35), constrained_layout=True)
    mesh = None
    for row, method in enumerate(methods):
        item = dataset[("main", method, -15)]
        left, right = axes[row]
        mesh = plot_image(left, item, FULL_BOUNDS, dynamic_range_db, show_zoom_box=True)
        plot_image(right, item, ZOOM_BOUNDS, dynamic_range_db)
        left.set_title(f"({chr(97 + 2 * row)}) {METHOD_LABELS[method]}: full FOV")
        right.set_title(f"({chr(98 + 2 * row)}) {METHOD_LABELS[method]}: zoom")
        left.set_ylabel("Range (m)")
        right.set_ylabel("Range (m)")
        if row == 2:
            left.set_xlabel("Angle (deg)")
            right.set_xlabel("Angle (deg)")
    assert mesh is not None
    add_shared_colorbar(fig, mesh, axes, dynamic_range_db)
    fig.text(
        0.5,
        1.002,
        f"Input SNR = −15 dB; FOV = 100°; display dynamic range = {dynamic_range_db:g} dB; + denotes ground truth",
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    save_figure(fig, output_dir / "Fig16_three_method_full_and_zoom")
    plt.close(fig)


def six_panel_figure(
    dataset: dict,
    output_dir: Path,
    dynamic_range_db: float,
    *,
    zoom: bool,
) -> None:
    bounds = ZOOM_BOUNDS if zoom else FULL_BOUNDS
    suffix = "zoom" if zoom else "full"
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.9), constrained_layout=True, sharex=True, sharey=True)
    mesh = None
    for index, (ax, method) in enumerate(zip(axes.ravel(), METHOD_ORDER)):
        item = dataset[("main", method, -15)]
        mesh = plot_image(ax, item, bounds, dynamic_range_db, show_zoom_box=False)
        ax.set_title(f"({chr(97 + index)}) {METHOD_LABELS[method]}")
        if index // 3 == 1:
            ax.set_xlabel("Angle (deg)")
        if index % 3 == 0:
            ax.set_ylabel("Range (m)")
    assert mesh is not None
    add_shared_colorbar(fig, mesh, axes, dynamic_range_db)
    descriptor = "zoom = 16° × 0.36 m" if zoom else "FOV = 100°"
    fig.text(
        0.5,
        1.002,
        f"Input SNR = −15 dB; {descriptor}; display dynamic range = {dynamic_range_db:g} dB; + denotes ground truth",
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    save_figure(fig, output_dir / f"Fig16_imaging_comparison_{suffix}")
    plt.close(fig)


def individual_figures(dataset: dict, output_dir: Path, dynamic_range_db: float) -> None:
    individual_dir = output_dir / "individual"
    for method in METHOD_ORDER:
        item = dataset[("main", method, -15)]
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), constrained_layout=True)
        mesh = plot_image(axes[0], item, FULL_BOUNDS, dynamic_range_db, show_zoom_box=True)
        plot_image(axes[1], item, ZOOM_BOUNDS, dynamic_range_db)
        axes[0].set_title(f"(a) {METHOD_LABELS[method]}: full FOV")
        axes[1].set_title(f"(b) {METHOD_LABELS[method]}: zoom (16° × 0.36 m)")
        for ax in axes:
            ax.set_xlabel("Angle (deg)")
            ax.set_ylabel("Range (m)")
        add_shared_colorbar(fig, mesh, axes, dynamic_range_db)
        fig.text(
            0.5,
            1.005,
            f"Input SNR = −15 dB; display dynamic range = {dynamic_range_db:g} dB; + denotes ground truth",
            ha="center",
            va="bottom",
            fontsize=8.0,
        )
        save_figure(fig, individual_dir / f"{method}_snr_m15_full_and_zoom")
        plt.close(fig)


def low_snr_figure(dataset: dict, output_dir: Path, dynamic_range_db: float) -> None:
    methods = ["dense", "cdm", "fdm"]
    snrs = [-10, -20, -30]
    fig, axes = plt.subplots(3, 3, figsize=(7.2, 7.0), constrained_layout=True, sharex=True, sharey=True)
    mesh = None
    for row, method in enumerate(methods):
        for col, snr in enumerate(snrs):
            item = dataset[("low_snr", method, snr)]
            ax = axes[row, col]
            mesh = plot_image(ax, item, FULL_BOUNDS, dynamic_range_db)
            ax.set_title(f"({chr(97 + row * 3 + col)}) {METHOD_LABELS[method]}, SNR = {snr} dB")
            if row == 2:
                ax.set_xlabel("Angle (deg)")
            if col == 0:
                ax.set_ylabel("Range (m)")
    assert mesh is not None
    add_shared_colorbar(fig, mesh, axes, dynamic_range_db)
    fig.text(
        0.5,
        1.002,
        f"FOV = 100°; display dynamic range = {dynamic_range_db:g} dB; + denotes ground-truth target",
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    save_figure(fig, output_dir / "Fig17_low_snr_comparison")
    plt.close(fig)


def write_figure_manifest(output_dir: Path, dynamic_range_db: float) -> None:
    rows = [
        {
            "figure": "Fig16_imaging_comparison_full_and_zoom",
            "panels": "6 methods x (full, zoom)",
            "snr_db": "-15",
            "fov_deg": "[-50, 50]",
            "range_m": "[4.80, 5.70]",
            "zoom": "angle [-8, 8] deg; range [4.96, 5.32] m; 16 deg x 0.36 m",
            "dynamic_range_db": dynamic_range_db,
        },
        {
            "figure": "Fig16_imaging_comparison_full",
            "panels": "6 methods",
            "snr_db": "-15",
            "fov_deg": "[-50, 50]",
            "range_m": "[4.80, 5.70]",
            "zoom": "none",
            "dynamic_range_db": dynamic_range_db,
        },
        {
            "figure": "Fig16_imaging_comparison_zoom",
            "panels": "6 methods",
            "snr_db": "-15",
            "fov_deg": "[-8, 8]",
            "range_m": "[4.96, 5.32]",
            "zoom": "16 deg x 0.36 m",
            "dynamic_range_db": dynamic_range_db,
        },
        {
            "figure": "Fig17_low_snr_comparison",
            "panels": "3 methods x 3 SNR values",
            "snr_db": "-10, -20, -30",
            "fov_deg": "[-50, 50]",
            "range_m": "[4.80, 5.70]",
            "zoom": "none",
            "dynamic_range_db": dynamic_range_db,
        },
    ]
    with (output_dir / "figure_manifest.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def verify_required_keys(dataset: dict) -> None:
    required = {("main", method, -15) for method in METHOD_ORDER}
    required |= {
        ("low_snr", method, snr)
        for method in ("dense", "cdm", "fdm")
        for snr in (-10, -20, -30)
    }
    missing = sorted(required - set(dataset))
    if missing:
        raise KeyError(f"Missing exported datasets: {missing}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    default_root = Path(__file__).resolve().parent
    parser.add_argument("--data-dir", type=Path, default=default_root / "data")
    parser.add_argument("--output-dir", type=Path, default=default_root / "figures")
    parser.add_argument("--dynamic-range-db", type=float, default=DEFAULT_DYNAMIC_RANGE_DB)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    apply_paper_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dataset = load_dataset(args.data_dir, args.dynamic_range_db)
    verify_required_keys(dataset)
    main_composite(dataset, args.output_dir, args.dynamic_range_db)
    three_method_composite(dataset, args.output_dir, args.dynamic_range_db)
    six_panel_figure(dataset, args.output_dir, args.dynamic_range_db, zoom=False)
    six_panel_figure(dataset, args.output_dir, args.dynamic_range_db, zoom=True)
    individual_figures(dataset, args.output_dir, args.dynamic_range_db)
    low_snr_figure(dataset, args.output_dir, args.dynamic_range_db)
    write_figure_manifest(args.output_dir, args.dynamic_range_db)
    print(f"Rendered reviewer-ready figures in {args.output_dir}")


if __name__ == "__main__":
    main()
