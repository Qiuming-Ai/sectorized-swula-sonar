from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
                     "font.size": 8, "savefig.dpi": 600})

fig, ax = plt.subplots(figsize=(7.15, 2.0))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
blue, orange, green, gray = "#DCEAF7", "#F8E1D2", "#DDEFE6", "#E8E8E8"
stages = [
    (0.015, 0.18, 0.205, 0.64, blue, "1  Spatial design",
     "Receive spacing and FOV\n→ grating directions\n→ six transmit sectors"),
    (0.265, 0.18, 0.205, 0.64, orange, "2  Waveform labels",
     "Simultaneous transmission\npacked FDM–LFM or CDM–ZC\nequal total pulse energy"),
    (0.515, 0.18, 0.205, 0.64, green, "3  Channel processing",
     "Sector matched filtering\nreceive DAS beamforming\nangle–range image"),
    (0.765, 0.18, 0.22, 0.64, gray, "4  Two checks",
     "Spatial: two-way alias masking\nWaveform: finite delay–Doppler\nand scene-density leakage"),
]
for x, y, w, h, color, title, body in stages:
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
                                linewidth=0.8, edgecolor="#555555", facecolor=color))
    ax.text(x + 0.014, y + h - 0.10, title, weight="bold", va="top")
    ax.text(x + 0.014, y + h - 0.24, body, va="top", linespacing=1.35)
for x in [0.225, 0.475, 0.725]:
    ax.add_patch(FancyArrowPatch((x, 0.50), (x + 0.035, 0.50), arrowstyle="-|>",
                                 mutation_scale=11, linewidth=1.0, color="#444444"))
ax.text(0.5, 0.06, "Transmit sectors address spatial aliasing; the waveform audit tests whether simultaneous sector echoes remain separable.",
        ha="center", fontsize=7, color="#444444")
fig.tight_layout(pad=0.15)
fig.savefig(FIGURES / "A18_workflow.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "A18_workflow.png", bbox_inches="tight")
print(FIGURES / "A18_workflow.pdf")
