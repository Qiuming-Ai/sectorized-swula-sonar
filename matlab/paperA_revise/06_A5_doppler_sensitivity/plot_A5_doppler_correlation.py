from __future__ import annotations

import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

with (RESULTS / "A5_correlation_profiles.csv").open(newline="", encoding="utf-8-sig") as f:
    profiles = list(csv.DictReader(f))
with (RESULTS / "A5_doppler_metrics.csv").open(newline="", encoding="utf-8-sig") as f:
    metrics = list(csv.DictReader(f))

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.7,
    "savefig.dpi": 600,
})

colors = {0.0: "#333333", 0.5: "#EE7733", 1.0: "#CC3311"}
styles = {0.0: "-", 0.5: "--", 1.0: "-."}
fig, axes = plt.subplots(1, 3, figsize=(7.15, 2.45))

for ax, method, title in zip(axes[:2], ["cdm", "fdm_current"], ["CDM-ZC", "FDM-LFM"]):
    for velocity in [0.0, 0.5, 1.0]:
        subset = [r for r in profiles if r["Method"] == method and abs(float(r["VelocityMps"])-velocity) < 1e-9]
        x = np.array([float(r["RangeOffsetM"]) for r in subset])
        y = np.maximum(-50, np.array([float(r["MagnitudeDb"]) for r in subset]))
        ax.plot(x, y, color=colors[velocity], linestyle=styles[velocity],
                label=f"v = {velocity:g} m/s")
    ax.set_xlim(-0.45, 0.45)
    ax.set_ylim(-50, 1)
    ax.set_xlabel("Apparent range offset (m)")
    ax.set_ylabel("Normalized correlation (dB)" if ax is axes[0] else "")
    ax.set_title(title, fontsize=8.5)
    ax.grid(True, color="#D9D9D9", linewidth=0.45)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].legend(frameon=False, loc="lower left")

for method, label, color, marker, linestyle in [
    ("cdm", "CDM-ZC", "#CC3311", "s", "--"),
    ("fdm_current", "FDM-LFM", "#0077BB", "^", "-."),
]:
    subset = sorted([r for r in metrics if r["Method"] == method], key=lambda r: float(r["VelocityMps"]))
    x = np.array([float(r["VelocityMps"]) for r in subset])
    y = np.array([float(r["MaxCrossAmbiguityDb"]) for r in subset])
    axes[2].plot(x, y, color=color, marker=marker, linestyle=linestyle, markersize=4, label=label)
axes[2].set_xlabel("Radial velocity (m/s)")
axes[2].set_ylabel("Worst cross-ambiguity peak (dB)")
axes[2].set_title("Inter-waveform leakage", fontsize=8.5)
axes[2].grid(True, color="#D9D9D9", linewidth=0.45)
axes[2].spines[["top", "right"]].set_visible(False)
axes[2].legend(frameon=False, loc="best")

for i, ax in enumerate(axes):
    ax.text(0.02, 0.96, f"({chr(97+i)})", transform=ax.transAxes, va="top", fontweight="bold")

fig.tight_layout(pad=0.7, w_pad=0.9)
fig.savefig(FIGURES / "A5_doppler_correlation.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "A5_doppler_correlation.png", bbox_inches="tight")
print(FIGURES / "A5_doppler_correlation.pdf")
