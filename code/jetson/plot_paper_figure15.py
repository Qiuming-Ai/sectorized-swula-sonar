"""Rebuild the manuscript's four-panel Jetson summary from released records.

This uses the manuscript layout (means and sample SD, plus explicit CUDA
allocation), not the three-panel distribution diagnostic. Vector SVG and
600-dpi PNG are exported; the archived manuscript PDF is not reproduced byte
for byte.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from verify_records import CODE, DEFAULT_RECORDS, METHODS, RELEASE, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--out", type=Path, default=RELEASE/"outputs"/"jetson"/"paper")
    args = parser.parse_args()
    report = verify(args.records)
    out = args.out
    outputs = [out/"figure15_jetson_resources.svg", out/"figure15_jetson_resources.png"]
    if any(p.exists() for p in outputs):
        raise FileExistsError("Choose a new output directory; existing figures are never overwritten")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.titlesize":10,
                         "svg.fonttype":"none","svg.hashsalt":"jetson-paper-figure15"})
    names = ["Dense","CDM","FDM"]
    colors = ["#345b84","#b85936","#328476"]
    fig, axes = plt.subplots(2,2,figsize=(6.4,4.6),layout="constrained")
    specs = (
        ("mean_ms","sample_SD_ms","(a) Complete reconstruction","Time (ms)"),
        ("explicit_CUDA_MiB",None,"(b) Explicit CUDA allocation","Memory (MiB)"),
        ("power_W_mean","power_W_sample_SD","(c) Module input power","Power (W)"),
        ("energy_J_per_block_mean","energy_J_per_block_sample_SD","(d) Energy per completed block","Energy (J/block)"),
    )
    for ax,(field,sd,title,ylabel) in zip(axes.flat,specs):
        values = [report["methods"][m][field] for m in METHODS]
        errors = [report["methods"][m][sd] if sd else 0. for m in METHODS]
        bars = ax.bar(names,values,color=colors,width=.58,yerr=errors,capsize=3)
        ax.set(title=title,ylabel=ylabel,ylim=(0,max(v+e for v,e in zip(values,errors))*1.22))
        ax.grid(axis="y",alpha=.18); ax.set_axisbelow(True)
        for bar,value in zip(bars,values):
            ax.annotate(f"{value:.3f}",(bar.get_x()+bar.get_width()/2,value),
                        xytext=(0,7),textcoords="offset points",ha="center",fontsize=9)
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(outputs[0],metadata={"Title":"Jetson reconstruction resource measurements","Date":None})
    fig.savefig(outputs[1],dpi=600,metadata={"Title":"Jetson reconstruction resource measurements"})
    plt.close(fig)
    print("Figure 15 layout rebuilt from verified timing, power, and memory records.")


if __name__ == "__main__":
    main()
