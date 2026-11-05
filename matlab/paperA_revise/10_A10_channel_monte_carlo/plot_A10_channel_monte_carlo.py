from __future__ import annotations

import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

with (RESULTS / "A10_summary_mean_ci.csv").open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

methods = ["dense", "cdm", "fdm_current"]
labels = {"dense": "Dense ULA", "cdm": "CDM SW-ULA", "fdm_current": "FDM SW-ULA"}
styles = {
    "dense": ("#333333", "o", "-"),
    "cdm": ("#CC3311", "s", "--"),
    "fdm_current": ("#0077BB", "^", "-."),
}
scenarios = [("point18", "18 point targets"), ("dense1000x50", "1000 scatterers / 50 range layers")]

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
    "savefig.dpi": 600,
})

fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.15), sharex="col")
for col, (scenario, title) in enumerate(scenarios):
    for method in methods:
        subset = sorted(
            [r for r in rows if r["Scenario"] == scenario and r["Method"] == method],
            key=lambda r: float(r["InputSNRDb"]),
        )
        x = np.array([float(r["InputSNRDb"]) for r in subset])
        color, marker, linestyle = styles[method]
        for row, mean_key, ci_key, ylabel in [
            (0, "MeanNMSEdB", "CI95NMSEdB", "Amplitude-aligned NMSE (dB)"),
            (1, "MeanCNR", "CI95CNR", "CNR"),
        ]:
            y = np.array([float(r[mean_key]) for r in subset])
            e = np.array([float(r[ci_key]) for r in subset])
            axes[row, col].errorbar(
                x, y, yerr=e, color=color, marker=marker, linestyle=linestyle,
                markersize=4, capsize=2.2, label=labels[method],
            )
            axes[row, col].set_ylabel(ylabel if col == 0 else "")
            axes[row, col].grid(True, color="#D9D9D9", linewidth=0.45)
            axes[row, col].spines[["top", "right"]].set_visible(False)
    axes[0, col].set_title(title, fontsize=8.5)
    axes[1, col].set_xlabel("Raw-channel input SNR (dB)")

axes[0, 0].legend(frameon=False, loc="best")
for i, ax in enumerate(axes.flat):
    ax.text(0.02, 0.96, f"({chr(97+i)})", transform=ax.transAxes, va="top", fontweight="bold")

fig.tight_layout(pad=0.8, w_pad=1.0, h_pad=0.8)
fig.savefig(FIGURES / "A10_channel_monte_carlo.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "A10_channel_monte_carlo.png", bbox_inches="tight")
print(FIGURES / "A10_channel_monte_carlo.pdf")
