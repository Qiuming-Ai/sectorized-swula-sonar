"""Plot already frozen density evidence, with no experiment rerun or tuning.

The main figure shows three independent scene values and their mean, not a
confidence interval. Full unbounded Student-t intervals stay in the supplement.
All generated outputs remain in this script's own integration directory.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import fitz

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "outputs/figures/figure12"
DATA = ROOT / "data/density"
RAW = DATA / "metrics.csv"
AGGREGATE = DATA / "aggregate_metrics.csv"
METHODS = [("dense", "Dense", "#0072B2", "o", "-"),
           ("cdm", "CDM", "#D55E00", "^", "--"),
           ("fdm", "FDM", "#009E73", "s", "-.")]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(RAW)
    a = pd.read_csv(AGGREGATE)
    selected = d[d.factor.isin(["scatter_count", "range_layers"])].copy()
    assert selected.metric_background_valid.eq(True).all()
    assert selected.support_leakage.between(0, 1).all()
    assert len(selected) == 99  # (six counts + five layer counts)*three scenes*three methods
    for _, q in selected.groupby(["factor", "factor_value"]):
        assert len(q) == 9 and q.scene_seed.nunique() == 3
        assert q.groupby("scene_seed").method.nunique().eq(3).all()
    means = selected.groupby(["factor", "factor_value", "method"], as_index=False).agg(
        rho_D0=("dense_magnitude_correlation", "mean"),
        outside_guard_power=("support_leakage", "mean"),
        available_background_min=("background_fraction", "min"),
        occupied_range_cells_min=("occupied_range_cells", "min"),
        occupied_range_cells_max=("occupied_range_cells", "max"))
    # Confirm the plotted means exactly reproduce the archived aggregation.
    matched = means.merge(a, on=["factor", "factor_value", "method"], validate="one_to_one")
    err_rho = float(np.max(abs(matched.rho_D0 - matched.dense_magnitude_correlation)))
    err_leak = float(np.max(abs(matched.outside_guard_power - matched.support_leakage)))
    assert err_rho < 1e-14 and err_leak < 1e-14
    means.to_csv(HERE / "density_main_values.csv", index=False)
    selected.to_csv(HERE / "density_main_individual_scenes.csv", index=False)
    plt.rcParams.update({"font.family": "Times New Roman", "font.size": 10.5,
        "axes.labelsize": 10.5, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 10.5, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": .7, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.75))
    settings = [("scatter_count", "Generated scatterer count (log scale)", "a", "b"),
                ("range_layers", "Requested range layers (log scale)", "c", "d")]
    for row, (factor, xlabel, letter_left, letter_right) in enumerate(settings):
        for col, (metric, label) in enumerate([
            ("dense_magnitude_correlation", "Correlation to noise-free Dense"),
            ("support_leakage", "Outside-guard power fraction")]):
            ax = axes[row, col]
            for method, name, colour, marker, style in METHODS:
                raw = selected[(selected.factor == factor) & (selected.method == method)]
                q = raw.groupby("factor_value", as_index=False)[metric].mean().sort_values("factor_value")
                ax.scatter(raw.factor_value, raw[metric], marker=marker, color=colour,
                           alpha=.32, s=15, linewidths=.35, zorder=2)
                ax.plot(q.factor_value, q[metric], color=colour, marker=marker, ls=style,
                        lw=1.25, ms=4.8, markerfacecolor="white", markeredgewidth=.85,
                        zorder=3, label=name)
            ax.set_xscale("log")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(label)
            ax.set_xlim((.85, 1200) if factor == "scatter_count" else (.85, 60))
            ticks = [1, 10, 100, 1000] if factor == "scatter_count" else [1, 5, 10, 20, 50]
            ax.set_xticks(ticks, [str(x) for x in ticks])
            ax.minorticks_off()
            ax.grid(axis="y", color="#BFBFBF", alpha=.38, lw=.5)
            ax.spines[["top", "right"]].set_visible(False)
            if col == 0:
                ax.set_ylim(.50, 1.04)
                ax.set_yticks([.5, .6, .7, .8, .9, 1.0])
            else:
                ax.set_ylim(0, .80)
                ax.set_yticks([0, .2, .4, .6, .8])
            ax.text(.015, 1.018, "(" + (letter_left if col == 0 else letter_right) + ")",
                    transform=ax.transAxes, va="bottom", fontsize=11)
    handles = [Line2D([0], [0], color=c, marker=m, ls=s, lw=1.25,
                     markerfacecolor="white", ms=5, label=n) for _, n, c, m, s in METHODS]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(.53, 1.00),
               frameon=False, handlelength=2.6, columnspacing=2.2)
    fig.subplots_adjust(left=.095, right=.985, top=.905, bottom=.115, wspace=.30, hspace=.43)
    outputs = []
    for ext in ["pdf", "svg", "png"]:
        path = HERE / ("figure_density_controls." + ext)
        fig.savefig(path, dpi=600)
        outputs.append(path)
    plt.close(fig)
    doc = fitz.open(outputs[0])
    spans = [s for b in doc[0].get_text("dict")["blocks"] if b.get("type") == 0
             for line in b["lines"] for s in line["spans"] if s.get("text", "").strip()]
    min_font = min(s["size"] for s in spans)
    doc[0].get_pixmap(matrix=fitz.Matrix(2, 2)).save(HERE / "density_main_render.png")
    receipt = {
        "source_metrics_csv": str(RAW), "source_sha256": sha(RAW),
        "source_aggregate_csv": str(AGGREGATE), "aggregate_sha256": sha(AGGREGATE),
        "rows": len(selected), "conditions": 11, "independent_generated_scenes_per_condition": 3,
        "added_noise": False, "methods_paired_within_scene": True,
        "display_statistics": "individual translucent scene points and arithmetic mean lines; no error bars",
        "supplement_statistics_unchanged": "unbounded Student-t 95% confidence intervals of each mean, n=3",
        "min_available_background_fraction": float(selected.background_fraction.min()),
        "all_background_metrics_valid": bool(selected.metric_background_valid.all()),
        "metric": "guard complement power / complete fixed ROI power; not fixed outer angular-band power",
        "guard": {"range_half_width_m": .06, "angle_half_width_deg": 1.0},
        "ROI": {"range_m": [30, 32.4], "angle_deg": [-50, 50]},
        "max_mean_disagreement": {"rho": err_rho, "outside_guard": err_leak},
        "source_font_min_pt": min_font,
        "safe_min_width_relative_to_6_5inch_canvas_for_8pt_font": float(8 / min_font),
        "raster_images_in_vector_pdf": len(doc[0].get_images()),
        "output_files": {str(p.name): sha(p) for p in outputs}}
    doc.close()
    (HERE / "DENSITY_MAIN_FIGURE_QA.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
