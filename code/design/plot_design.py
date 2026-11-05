"""Portable design-figure redraws from documented legacy inputs.

The outputs preserve the stated numerical or diagram content, not the
original PDF byte streams. These are not the final implemented-pulse or
complete-bank scene experiments.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as VectorPath
from matplotlib.patches import PathPatch
import numpy as np
from scipy.signal import fftconvolve


ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "data/design"
SPATIAL = ROOT / "data/spatial"
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
STYLES = ["-", "--", "-.", ":", (0, (5, 2)), (0, (3, 1, 1, 1))]
plt.rcParams.update({"font.family": "serif", "font.size": 9, "pdf.fonttype": 42})


def style_axes(ax, xlabel, ylabel, xlim, ylim):
    ax.set(xlabel=xlabel, ylabel=ylabel, xlim=xlim, ylim=ylim)
    ax.grid(alpha=.2)


def figure2():
    a = np.loadtxt(DESIGN / "figure2_legacy_psf_slices.csv", delimiter=",", skiprows=1)
    assert a.shape == (1701, 8)
    with np.errstate(divide="ignore"):
        db = 20 * np.log10(a[:, 1:])
    fig, axes = plt.subplots(2, 1, figsize=(5.7, 4.4), sharex=True, layout="constrained")
    for k, label in enumerate(["Rx focus 0°", "Rx focus −6°", "Rx focus +6°"]):
        axes[0].plot(a[:, 0], db[:, k], label=label, color=PALETTE[k], ls=STYLES[k])
        axes[1].plot(a[:, 0], db[:, 4 + k], label=label, color=PALETTE[k], ls=STYLES[k])
    axes[0].plot(a[:, 0], db[:, 3], label="Tx focus 0°", color=PALETTE[3], ls=STYLES[3])
    for ax, title in zip(axes, ("(a) Transmit and receive slices", "(b) Two-way product slices")):
        style_axes(ax, "Target bearing (deg)" if ax is axes[1] else "",
                   "Normalized magnitude (dB)", (-40, 40), (-60, 1))
        ax.set_title(title)
        ax.legend(loc="lower center", ncol=2, fontsize=8)
    return fig


def figure3():
    a = np.loadtxt(SPATIAL / "Figure03_transmit_beampattern_data.csv", delimiter=",", skiprows=1)
    cfg = json.loads((DESIGN / "A1_parameters.json").read_text(encoding="utf-8"))
    assert a.shape == (2001, 7)
    fig, ax = plt.subplots(figsize=(3.35, 2.6), layout="constrained")
    for k, focus in enumerate(cfg["focus_angles_deg"]):
        ax.plot(a[:, 0], a[:, k + 1], color=PALETTE[k], ls=STYLES[k], label=f"{focus:g}°")
    for edge in np.unique(cfg["sector_edges_deg"]):
        ax.axvline(edge, color=".6", lw=.55, ls=":")
    style_axes(ax, "Bearing (deg)", "Normalized magnitude (dB)", (-50, 50), (-45, 2))
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.25), ncol=3, frameon=False, fontsize=8)
    return fig


def legacy_fdm_design():
    """Original zero-Doppler prototype, including the historical cosine envelope."""
    duration, bandwidth, fs = .5e-3, 20e3, 2e6
    t = np.arange(round(duration * fs)) / fs
    w = [np.ones(t.size), np.hanning(t.size), np.hamming(t.size),
         np.blackman(t.size), np.kaiser(t.size, 6)]
    df_ratio = np.linspace(0, 2.5, 1001)
    bank = np.asarray(w) * np.cos(2 * np.pi * t / duration) * np.exp(
        1j * np.pi * (bandwidth / duration) * t**2,
    )
    offsets = np.exp(2j * np.pi * bandwidth * df_ratio[:, None] * t)
    ratio_db = []
    for reference in bank:
        received = offsets * reference
        cross = fftconvolve(received, reference.conj()[None, ::-1], axes=1)
        auto_peak = np.vdot(reference, reference).real
        ratio_db.append(20 * np.log10(np.maximum(abs(cross).max(axis=1) / auto_peak, 1e-5)))
    result = np.asarray(ratio_db)
    # Frozen spot checks at df/B = 2.0, rounded to the historical chart's
    # readable precision. These do not validate the implemented waveform bank.
    np.testing.assert_allclose(result[[0, 2, 4], 800], [-24.5, -42.6, -56.8], atol=.07)
    return df_ratio, result


def legacy_zc_design():
    """Ideal periodic upper bound used in the historical CDM example."""
    length = np.unique(np.rint(np.logspace(np.log10(16), np.log10(8192), 100))).astype(int)
    beam_count = np.array([2, 4, 6, 8])
    pulse_ms = length / 100.0  # L / 100 kchip/s in ms
    ratios = 20 * np.log10(beam_count[:, None] - 1) - 10 * np.log10(length)
    return pulse_ms, ratios, beam_count


def figure4a():
    df_ratio, result = legacy_fdm_design()
    fig, ax = plt.subplots(figsize=(4.25, 3), layout="constrained")
    for k, label in enumerate(["Rect", "Hann", "Hamming", "Blackman", "Kaiser β=6"]):
        ax.plot(df_ratio, result[k], color=PALETTE[k], ls=STYLES[k], label=label)
    style_axes(ax, "Subband spacing / bandwidth", "Peak amplitude ratio (dB)", (0, 2), (-60, 2))
    ax.set_title("FDM prototype, T = 0.5 ms, B = 20 kHz")
    ax.legend(fontsize=8)
    return fig


def figure4b():
    pulse_ms, ratios, beam_count = legacy_zc_design()
    fig, ax = plt.subplots(figsize=(4.25, 3), layout="constrained")
    for k, m in enumerate(beam_count):
        ax.plot(pulse_ms, ratios[k], color=PALETTE[k], ls=STYLES[k], label=f"M = {m}")
    style_axes(ax, "Pulse duration (ms)", "Peak amplitude ratio (dB)", (0, 10), (-42, 6))
    ax.set_title("Ideal periodic CDM-ZC bound, B = 100 kHz")
    ax.legend(fontsize=8)
    return fig


def figure7():
    d = json.loads((DESIGN / "figure7_pipeline_geometry.json").read_text(encoding="utf-8"))
    assert len(d["paths"]) == 22 and len(d["nodes"]) == 8
    w, h = d["width_pt"], d["height_pt"]
    fig = plt.figure(figsize=(w / 72, h / 72))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, w), ylim=(0, h))
    ax.axis("off")
    for path in d["paths"]:
        points, codes = [path["items"][0][1]], [VectorPath.MOVETO]
        for primitive, _start, *end in path["items"]:
            points.extend(end)
            codes.extend([VectorPath.LINETO] if primitive == "l" else [VectorPath.CURVE4] * 3)
        if path["closed"]:
            points.append(points[0])
            codes.append(VectorPath.CLOSEPOLY)
        ax.add_patch(PathPatch(
            VectorPath(points, codes), facecolor=path["fill_rgb"] or "none",
            edgecolor=path["stroke_rgb"] or "none", lw=path["width_pt"],
        ))
    for node in d["nodes"]:
        x, y = node["center_pt"]
        ax.text(x, y, "\n".join(node["lines"]), ha="center", va="center", fontsize=9)
    return fig


def figure8():
    with (SPATIAL / "A14_angular_profiles.csv").open(encoding="utf-8-sig", newline="") as f:
        rows = [r for r in csv.DictReader(f) if float(r["TargetAngleDeg"]) == 8.]
    selections = [
        [("pure_receive", "dense", "Dense ULA", "k", "-"),
         ("pure_receive", "sw_ula", "SW-ULA (common Rx)", PALETTE[0], "--")],
        [("effective_angular_psf", "dense", "Dense ULA", "k", "-"),
         ("effective_angular_psf", "cdm", "CDM effective PSF", PALETTE[1], "--"),
         ("effective_angular_psf", "fdm_current", "Packed FDM effective PSF", PALETTE[0], "-.")],
    ]
    fig, axes = plt.subplots(2, 1, figsize=(3.35, 4.55), layout="constrained")
    for ax, selection, label in zip(axes, selections, ("(a) Receive only", "(b) Effective PSF")):
        for kind, method, series_label, color, ls in selection:
            xy = sorted(
                (float(r["ScanAngleDeg"]), max(-55, float(r["MagnitudeDb"])))
                for r in rows if r["ResponseType"] == kind and r["Method"] == method
            )
            assert len(xy) == 6001
            a = np.asarray(xy)
            ax.plot(a[:, 0], a[:, 1], color=color, ls=ls, label=series_label)
        style_axes(ax, "Scan bearing (deg)", "Magnitude (dB)", (-50, 50), (-55, 2))
        ax.set_title(label)
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.27), frameon=False, fontsize=8)
    return fig


FIGURES = {"2": figure2, "3": figure3, "4a": figure4a, "4b": figure4b,
           "7": figure7, "8": figure8}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures", nargs="+", choices=FIGURES, default=list(FIGURES))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify-only", action="store_true",
                        help="Build each vector PDF in memory without writing outputs.")
    args = parser.parse_args()
    if not args.verify_only and args.output is None:
        parser.error("provide --output or --verify-only")
    if not args.verify_only:
        args.output.mkdir(parents=True, exist_ok=True)
    for key in args.figures:
        fig = FIGURES[key]()
        if args.verify_only:
            stream = io.BytesIO()
            fig.savefig(stream, format="pdf")
            assert stream.getvalue().startswith(b"%PDF")
        else:
            for suffix in (".pdf", ".png"):
                fig.savefig(args.output / f"figure{key}_design_redraw{suffix}", dpi=300)
        plt.close(fig)
        print(f"Figure {key}: rendered")


if __name__ == "__main__":
    main()
