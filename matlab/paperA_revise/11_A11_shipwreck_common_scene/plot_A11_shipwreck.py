from __future__ import annotations

from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[3]
parser = argparse.ArgumentParser(description="Plot historical MATLAB A11 outputs, not v08 Figure 13.")
parser.add_argument("--results", type=Path, default=ROOT / "outputs/matlab/A11_paper_archive")
RESULTS = parser.parse_args().results
FIGURES = RESULTS / "figures"

d = loadmat(RESULTS / "A11_plot_data.mat", squeeze_me=True)
FIGURES.mkdir(exist_ok=True)
angle = np.asarray(d["angle_deg"]).squeeze()
range_m = np.asarray(d["range_m"]).squeeze()
truth = np.asarray(d["truth_magnitude"], dtype=float)
images = np.asarray(d["image_magnitude"], dtype=float)

panels = [truth, images[:, :, 0], images[:, :, 1], images[:, :, 2]]
labels = [
    "Ray-traced reference (display enhanced)",
    "Dense ULA",
    "CDM SW-ULA",
    "FDM SW-ULA",
]

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "savefig.dpi": 600,
})

fig, axes = plt.subplots(2, 2, figsize=(7.15, 5.45), sharex=True, sharey=True)
last = None
for idx, (ax, panel, label) in enumerate(zip(axes.flat, panels, labels)):
    panel_db = 20 * np.log10(panel / (np.max(panel) + np.finfo(float).eps) + np.finfo(float).eps)
    last = ax.imshow(
        panel_db,
        extent=[angle.min(), angle.max(), range_m.min(), range_m.max()],
        origin="lower",
        aspect="auto",
        vmin=-40,
        vmax=0,
        cmap="cividis",
        interpolation="nearest",
    )
    ax.text(0.02, 0.96, f"({chr(97+idx)}) {label}", transform=ax.transAxes,
            va="top", ha="left", color="white", fontweight="bold",
            bbox=dict(facecolor="black", alpha=0.42, edgecolor="none", pad=1.5))
    ax.set_xlim(-50, 50)
    ax.set_ylim(range_m.min(), range_m.max())
    if idx // 2 == 1:
        ax.set_xlabel("Receive-array angle, asin(x/r) (deg)")
    if idx % 2 == 0:
        ax.set_ylabel("Slant range (m)")

fig.subplots_adjust(left=0.09, right=0.88, bottom=0.10, top=0.98, wspace=0.08, hspace=0.10)
cax = fig.add_axes([0.90, 0.16, 0.018, 0.70])
cbar = fig.colorbar(last, cax=cax)
cbar.set_label("Normalized display level (dB)")
fig.savefig(FIGURES / "A11_shipwreck_common_scene.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "A11_shipwreck_common_scene.png", bbox_inches="tight")
print(FIGURES / "A11_shipwreck_common_scene.pdf")

# A receive-array polar-view companion is exported for visual comparison with
# sector_plot in the original Exp-8 script. With sin(theta)=x/r, the displayed
# lateral coordinate equals x, whereas r*cos(theta)=sqrt(y^2+z^2) is an
# unresolved-plane radial coordinate rather than the physical forward z axis.
theta, radial = np.meshgrid(np.deg2rad(angle), range_m)
x_sector = radial * np.sin(theta)
z_sector = radial * np.cos(theta)
sector_fig, sector_axes = plt.subplots(2, 2, figsize=(7.15, 3.95), sharex=True, sharey=True)
sector_last = None
for idx, (ax, panel, label) in enumerate(zip(sector_axes.flat, panels, labels)):
    panel_db = 20 * np.log10(panel / (np.max(panel) + np.finfo(float).eps) + np.finfo(float).eps)
    sector_last = ax.pcolormesh(
        x_sector, z_sector, panel_db, shading="auto", vmin=-40, vmax=0,
        cmap="cividis", rasterized=True,
    )
    ax.text(0.02, 0.96, f"({chr(97+idx)}) {label}", transform=ax.transAxes,
            va="top", ha="left", color="white", fontweight="bold",
            bbox=dict(facecolor="black", alpha=0.42, edgecolor="none", pad=1.5))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-95, 95)
    ax.set_ylim(26, 122)
    if idx // 2 == 1:
        ax.set_xlabel("Lateral coordinate x (m)")
    if idx % 2 == 0:
        ax.set_ylabel(r"Radial coordinate $R_\perp$ (m)")
sector_fig.subplots_adjust(left=0.085, right=0.88, bottom=0.12, top=0.98,
                           wspace=0.08, hspace=0.16)
sector_cax = sector_fig.add_axes([0.90, 0.17, 0.018, 0.68])
sector_bar = sector_fig.colorbar(sector_last, cax=sector_cax)
sector_bar.set_label("Normalized display level (dB)")
sector_fig.savefig(FIGURES / "A11_shipwreck_common_scene_sector.pdf", bbox_inches="tight")
sector_fig.savefig(FIGURES / "A11_shipwreck_common_scene_sector.png", bbox_inches="tight")
print(FIGURES / "A11_shipwreck_common_scene_sector.pdf")

# Geometry audit: verify that the scene retains the legacy pose and the
# 45-deg-elevation line of sight before interpreting the range--azimuth image.
ship = np.genfromtxt(RESULTS / "A11_ship_facets.csv", delimiter=",", names=True)
seabed = np.genfromtxt(RESULTS / "A11_seabed_scatterers.csv", delimiter=",", names=True)
ship_step = max(1, len(ship) // 18000)
seabed_step = max(1, len(seabed) // 12000)
ship_q = ship[::ship_step]
seabed_q = seabed[::seabed_step]

audit = plt.figure(figsize=(7.15, 3.15))
ax3 = audit.add_subplot(1, 2, 1, projection="3d")
ax3.scatter(seabed_q["X"], seabed_q["Z"], seabed_q["Y"], s=0.12,
            c="#9aa0a6", alpha=0.22, rasterized=True, label="Legacy seabed")
ax3.scatter(ship_q["X"], ship_q["Z"], ship_q["Y"], s=0.25,
            c=ship_q["RangeM"], cmap="viridis", alpha=0.65, rasterized=True,
            label="ship.stl facets")
ax3.scatter([0], [0], [0], marker="^", s=28, c="#d95f02", label="Sonar origin")
look = 42 / np.sqrt(2)
ax3.quiver(0, 0, 0, 0, look, look, color="#d95f02", linewidth=1.2,
           arrow_length_ratio=0.08)
ax3.set_xlabel("Lateral x (m)")
ax3.set_ylabel("Forward z (m)")
ax3.set_zlabel("Vertical y (m)")
ax3.view_init(elev=22, azim=-58)
ax3.set_title("(a) Legacy 3-D pose and inactive 45° preview direction", loc="left", fontsize=8)
ax3.legend(loc="upper left", fontsize=6, frameon=False)

ax2 = audit.add_subplot(1, 2, 2)
ax2.scatter(seabed_q["AngleDeg"], seabed_q["RangeM"], s=0.10,
            c="#9aa0a6", alpha=0.22, rasterized=True, label="Legacy seabed")
points = ax2.scatter(ship_q["AngleDeg"], ship_q["RangeM"], s=0.22,
                     c=ship_q["ElevationDeg"], cmap="coolwarm",
                     vmin=-35, vmax=35, alpha=0.68, rasterized=True,
                     label="Visible ship facets")
ax2.set_xlim(-50, 50)
ax2.set_ylim(range_m.min(), range_m.max())
ax2.set_xlabel("Receive-array direction-cosine angle, asin(x/r) (deg)")
ax2.set_ylabel("Slant range (m)")
ax2.set_title("(b) Common 2-D projection used by all methods", loc="left", fontsize=8)
ax2.grid(color="0.88", linewidth=0.45)
cbar2 = audit.colorbar(points, ax=ax2, pad=0.02, fraction=0.05)
cbar2.set_label("Unresolved elevation (deg)")
audit.subplots_adjust(left=0.06, right=0.97, bottom=0.16, top=0.94, wspace=0.28)
audit.savefig(FIGURES / "A11_geometry_audit.pdf", bbox_inches="tight")
audit.savefig(FIGURES / "A11_geometry_audit.png", bbox_inches="tight")
print(FIGURES / "A11_geometry_audit.pdf")
