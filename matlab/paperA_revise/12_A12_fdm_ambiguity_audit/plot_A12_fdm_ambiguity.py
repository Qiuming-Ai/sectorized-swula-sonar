from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)


def read_csv(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


summary = read_csv("A12_doppler_summary.csv")
profiles = read_csv("A12_delay_profiles.csv")
spectra = read_csv("A12_spectra.csv")

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.8,
        "axes.linewidth": 0.7,
        "savefig.dpi": 600,
    }
)

blue = "#0072B2"
orange = "#D55E00"
green = "#009E73"
purple = "#CC79A7"
gray = "#555555"

fig, axes = plt.subplots(2, 2, figsize=(7.15, 4.75))

# (a) Full finite-pulse delay profile for the worst FDM pair at v=0.
fdm = [
    row
    for row in profiles
    if row["Method"] == "fdm_current" and abs(float(row["VelocityMps"])) < 1e-12
]
x = np.asarray([float(row["DelayMs"]) for row in fdm])
y = np.maximum(-100, np.asarray([float(row["MagnitudeDb"]) for row in fdm]))
axes[0, 0].plot(x, y, color=blue, linewidth=1.1)
axes[0, 0].axvspan(-0.6, 0.6, color="#D9EAF7", alpha=0.8, label="evaluated image window")
peak_idx = int(np.argmax(y))
axes[0, 0].plot(x[peak_idx], y[peak_idx], "o", color=orange, markersize=4)
axes[0, 0].annotate(
    f"{y[peak_idx]:.1f} dB at {x[peak_idx]:.2f} ms",
    xy=(x[peak_idx], y[peak_idx]),
    xytext=(-3.4, -14),
    arrowprops={"arrowstyle": "->", "color": gray, "lw": 0.7},
    fontsize=7,
)
axes[0, 0].set(xlabel="Relative delay (ms)", ylabel="Normalized magnitude (dB)", ylim=(-100, 2))
axes[0, 0].set_title("Complete finite-pulse delay search")

# (b) Same profile restricted to the range span used by the point-scene audit.
local = np.abs(x) <= 0.6 + 1e-12
axes[0, 1].plot(x[local], y[local], color=blue, linewidth=1.1)
local_idx = np.flatnonzero(local)[int(np.argmax(y[local]))]
axes[0, 1].plot(x[local_idx], y[local_idx], "o", color=orange, markersize=4)
axes[0, 1].annotate(
    f"{y[local_idx]:.1f} dB",
    xy=(x[local_idx], y[local_idx]),
    xytext=(-0.42, -61),
    arrowprops={"arrowstyle": "->", "color": gray, "lw": 0.7},
    fontsize=7,
)
axes[0, 1].set(xlabel="Relative delay (ms)", ylabel="Normalized magnitude (dB)", ylim=(-100, -55))
axes[0, 1].set_title("Restricted ±0.6-ms processing window")

# (c) Doppler curves: complete search versus restricted window.
style = {
    ("fdm_current", "full"): (blue, "o", "-", "FDM, complete delay"),
    ("fdm_current", "local"): (blue, "^", "--", "FDM, ±0.6 ms"),
    ("cdm", "full"): (orange, "s", "-", "CDM, complete delay"),
    ("cdm", "local"): (orange, "v", "--", "CDM, ±0.6 ms"),
}
for method in ["fdm_current", "cdm"]:
    block = sorted(
        [row for row in summary if row["Method"] == method],
        key=lambda row: float(row["VelocityMps"]),
    )
    velocity = np.asarray([float(row["VelocityMps"]) for row in block])
    for scope, field in [("full", "WorstFullPeakDb"), ("local", "WorstProcessingWindowPeakDb")]:
        color, marker, linestyle, label = style[(method, scope)]
        values = np.asarray([float(row[field]) for row in block])
        axes[1, 0].plot(
            velocity,
            values,
            color=color,
            marker=marker,
            linestyle=linestyle,
            linewidth=1,
            markersize=3.5,
            label=label,
        )
axes[1, 0].set(xlabel="Radial velocity (m/s)", ylabel="Worst cross-ambiguity peak (dB)")
axes[1, 0].set_title("Search-domain dependence of reported leakage")
axes[1, 0].legend(frameon=False, ncol=2, loc="lower center")

# (d) Implemented FDM spectra and RF subband edges.
for waveform, color in zip(range(1, 7), [blue, orange, green, purple, "#E69F00", "#56B4E9"]):
    block = [
        row
        for row in spectra
        if row["Method"] == "fdm_current" and int(float(row["Waveform"])) == waveform
    ]
    frequency = np.asarray([float(row["RFFrequencyKHz"]) for row in block])
    magnitude = np.asarray([float(row["MagnitudeDb"]) for row in block])
    sector = float(block[0]["SectorAngleDeg"])
    axes[1, 1].plot(frequency, magnitude, color=color, linewidth=0.85, label=f"{sector:g}°")
axes[1, 1].set(xlabel="RF frequency (kHz)", ylabel="Normalized spectrum (dB)", xlim=(395, 505), ylim=(-70, 2))
axes[1, 1].set_title("Packed 20-kHz FDM–LFM waveforms")
axes[1, 1].legend(frameon=True, facecolor="white", edgecolor="none", ncol=3,
                  loc="lower center", bbox_to_anchor=(0.5, -0.02), title="Tx sector")

for index, ax in enumerate(axes.flat):
    ax.grid(True, color="#D9D9D9", linewidth=0.45)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.02, 0.96, f"({chr(97 + index)})", transform=ax.transAxes, va="top", fontweight="bold")

fig.tight_layout(pad=0.8, w_pad=1.0, h_pad=1.0)
fig.savefig(FIGURES / "A12_fdm_ambiguity_audit.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "A12_fdm_ambiguity_audit.png", bbox_inches="tight")
print(FIGURES / "A12_fdm_ambiguity_audit.pdf")
