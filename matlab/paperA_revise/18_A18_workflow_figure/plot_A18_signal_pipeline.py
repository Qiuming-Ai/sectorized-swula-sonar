from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 7.5,
})

labels = [
    "Signal\ntransmission",
    "Channel\npropagation",
    "Matched\nfiltering",
    "Shift and\npadding",
    "FFT",
    "Phase\ncompensation",
    "Weighted\nsum",
    "IFFT",
]
widths = [0.105, 0.115, 0.105, 0.105, 0.070, 0.125, 0.095, 0.070]
colors = ["#DCEAF7", "#F2E6D8", "#DDEFE6", "#E8E8E8",
          "#DCEAF7", "#F2E6D8", "#DDEFE6", "#E8E8E8"]

fig, ax = plt.subplots(figsize=(7.15, 0.78))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

x = 0.015
gap = 0.021
for index, (label, width, color) in enumerate(zip(labels, widths, colors)):
    ax.add_patch(FancyBboxPatch(
        (x, 0.22), width, 0.56,
        boxstyle="round,pad=0.007,rounding_size=0.012",
        facecolor=color, edgecolor="#4D4D4D", linewidth=0.75,
    ))
    ax.text(x + width / 2, 0.50, label, ha="center", va="center", linespacing=1.1)
    if index < len(labels) - 1:
        ax.add_patch(FancyArrowPatch(
            (x + width + 0.003, 0.50), (x + width + gap - 0.002, 0.50),
            arrowstyle="-|>", mutation_scale=9, linewidth=0.85, color="#444444",
        ))
    x += width + gap

fig.tight_layout(pad=0.05)
fig.savefig(FIGURES / "A18_signal_pipeline.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "A18_signal_pipeline.png", dpi=600, bbox_inches="tight")
print(FIGURES / "A18_signal_pipeline.pdf")
