"""Render completed, valid Jetson power windows without smoothing or imputation.

Acquisition sources, engineering engines, receipts and frozen manifests are
read-only. The independent analysis code and input hashes accompany each export.
Requires all 3 sessions x 3 methods before creating any result figure.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics
import sys
from pathlib import Path

METHODS = ("dense", "cdm", "fdm")
SESSIONS = (0, 1, 2)
LABELS = {"dense": "Dense LFM", "cdm": "CDM-ZC", "fdm": "FDM-LFM"}
# Okabe-Ito: every method also has a distinct shape and line style.
COLORS = {"dense": "#0072B2", "cdm": "#D55E00", "fdm": "#009E73"}
MARKERS = {"dense": "o", "cdm": "^", "fdm": "s"}
LINESTYLES = {"dense": "-", "cdm": "--", "fdm": ":"}
PHASE_COLORS = {"Idle": "#F0F0F0", "Warm-up": "#FFF1D6", "Measured window": "#E2EFEB"}
EXPECTED_RELATIVE_DURATION = 120.0
OUTPUT_STEMS = ("jetson_summary", "module_power_energy", "power_time_all_windows")
RELEASE_ROOT = Path(__file__).resolve().parents[3]


def read_json(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Complete nine-window campaign required; missing {path}")
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1_048_576), b""):
            h.update(block)
    return h.hexdigest()


def finite_number(value, label, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Nonfinite/missing {label}")
    if positive and value <= 0:
        raise ValueError(f"Nonpositive {label}")
    return float(value)


def exact_integer(value, label, positive=False):
    if isinstance(value, bool) or not isinstance(value, int) or (positive and value <= 0):
        raise ValueError(f"Invalid integer {label}")
    return value


def close(a, b, label, tolerance=1e-10):
    if not math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError(f"Inconsistent {label}: {a} != {b}")


def one_event(invocation, name):
    selected = [e for e in invocation.get("events", []) if e.get("event") == name]
    if len(selected) != 1:
        raise ValueError(f"Exactly one {name} runner event required")
    return selected[0]


def load_window(run_root, session, method):
    directory = run_root / "power" / f"session{session}" / method
    paths = (
        directory / "POWER_INTEGRAL.json",
        directory / "power_samples.json",
        directory / "runner" / "window.json",
        directory / "invocation.json",
    )
    integral, samples, window, invocation = [read_json(p) for p in paths]
    if integral.get("valid") is not True or integral.get("reasons"):
        raise ValueError(f"Invalid integrated window, no result figure generated: {directory}")
    if window.get("workload_valid") is not True or invocation.get("returncode") != 0:
        raise ValueError(f"Invalid/incomplete runner window: {directory}")
    if invocation.get("whole_runner_timeout"):
        raise ValueError(f"Runner timeout retained, plotting rejected: {directory}")
    if integral.get("idle_subtracted") is not False:
        raise ValueError("Primary graph must use module-total power without idle subtraction")
    if window.get("method") != method or window.get("session") != session:
        raise ValueError(f"Window identity mismatch: {directory}")
    t0 = exact_integer(window.get("t0_ns"), "t0_ns", True)
    t1 = exact_integer(window.get("t1_ns"), "t1_ns", True)
    count = exact_integer(window.get("completed_blocks"), "completed_blocks", True)
    if t1 <= t0 or integral.get("t0_ns") != t0 or integral.get("t1_ns") != t1 or integral.get("completed_blocks") != count:
        raise ValueError(f"Window/integral endpoint or completed-count mismatch: {directory}")
    duration = (t1 - t0) / 1e9
    target = window.get("target_submission_ns")
    if target != 120_000_000_000 or duration < EXPECTED_RELATIVE_DURATION:
        raise ValueError(f"Unexpected submission/drain window duration: {directory}")
    average = finite_number(integral.get("average_module_power_W"), "average module power", True)
    energy = finite_number(integral.get("energy_J"), "window energy", True)
    per_block = finite_number(integral.get("energy_J_per_block"), "energy per block", True)
    close(integral.get("actual_window_seconds"), duration, "actual duration")
    close(average * duration, energy, "average power * actual duration")
    close(per_block * count, energy, "energy per block * completed count")

    idle = one_event(invocation, "idle_begin")
    warm = one_event(invocation, "warmup_begin")
    start = one_event(invocation, "window_start")
    end = one_event(invocation, "window_end")
    idle_ns = exact_integer(idle.get("monotonic_ns"), "idle event", True)
    warm_ns = exact_integer(warm.get("monotonic_ns"), "warmup event", True)
    if not idle_ns < warm_ns < t0 < t1:
        raise ValueError("Runner phase timestamps do not have the declared order")
    if start.get("t0_ns") != t0 or start.get("monotonic_ns") != t0 or end.get("t0_ns") != t0 or end.get("t1_ns") != t1 or end.get("monotonic_ns") != t1:
        raise ValueError("Runner phase events and window timestamps disagree")
    if end.get("completed_blocks") != count or end.get("workload_valid") is not True:
        raise ValueError("Runner window_end count/status mismatch")

    rows = samples.get("rows", [])
    if not rows:
        raise ValueError("No raw telemetry samples")
    times = [exact_integer(r.get("timestamp_ns"), "telemetry timestamp", True) for r in rows]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("Telemetry must retain its strictly increasing original order")
    power_rows = [r for r in rows if r.get("power_W") is not None]
    if not power_rows:
        raise ValueError("No valid current VDD_IN reads")
    for r in power_rows:
        finite_number(r["power_W"], "raw power")
        if r["power_W"] < 0:
            raise ValueError("Negative current VDD_IN read")
    before = [r for r in power_rows if r["timestamp_ns"] <= t0]
    after = [r for r in power_rows if r["timestamp_ns"] >= t1]
    if not before or not after:
        raise ValueError("Unbracketed raw power endpoints")
    region = [r for r in rows if before[-1]["timestamp_ns"] <= r["timestamp_ns"] <= after[0]["timestamp_ns"]]
    if any(r.get("power_W") is None for r in region):
        raise ValueError("Missing raw power reading in integrated window")
    if any(b["timestamp_ns"] - a["timestamp_ns"] > 500_000_000 for a, b in zip(region, region[1:])):
        raise ValueError("Gap over 500 ms in integrated window")
    temperature_rows = [r for r in rows if r.get("max_temperature_C") is not None]
    for r in temperature_rows:
        finite_number(r["max_temperature_C"], "maximum available temperature")
    if any(r["max_temperature_C"] >= 80 for r in region if r.get("max_temperature_C") is not None):
        raise ValueError("Safety-stop temperature in measured window")
    return {
        "session": session, "method": method, "directory": directory,
        "t0_ns": t0, "t1_ns": t1, "idle_ns": idle_ns, "warm_ns": warm_ns,
        "duration_s": duration, "completed_blocks": count,
        "average_power_W": average, "energy_per_block_J": per_block, "energy_J": energy,
        "rows": rows, "power_rows": power_rows, "temperature_rows": temperature_rows,
        "auxiliary_reader_errors": samples.get("reader_errors", []),
        "sources": paths,
    }


def load_campaign(run_root):
    windows = [load_window(run_root, s, m) for s in SESSIONS for m in METHODS]
    if len(windows) != 9:
        raise ValueError("Nine complete independent windows are required")
    return windows


def load_latency_campaign(run_root):
    """Keep all 600 measured calls per method; warmup exclusion is by protocol."""
    values = {m: [] for m in METHODS}
    sessions = {m: [] for m in METHODS}
    sources = []
    for session in SESSIONS:
        for method in METHODS:
            path = run_root / "latency" / f"session{session}" / method / "runner" / "latency_raw.csv"
            if not path.is_file():
                raise FileNotFoundError(f"Complete 1800-call latency catalog required; missing {path}")
            with path.open("r", newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            if len(rows) != 250 or any(r.get("phase") not in ("warmup", "measure") for r in rows):
                raise ValueError(f"Expected exactly 50 warmup + 200 measure rows: {path}")
            warm = [r for r in rows if r["phase"] == "warmup"]
            measured = [r for r in rows if r["phase"] == "measure"]
            if len(warm) != 50 or len(measured) != 200:
                raise ValueError(f"Latency phase count mismatch: {path}")
            for phase, expected in ((warm, 50), (measured, 200)):
                if [int(r["block_id"]) for r in phase] != list(range(expected)):
                    raise ValueError(f"Original block IDs must be complete and ordered: {path}")
                for r in phase:
                    if int(r["session"]) != session or r["method"] != method or r["status"] != "ok":
                        raise ValueError(f"Invalid/mislabeled latency call retained, graph rejected: {path}")
                    if r["shape_pass"] != "1" or r["checkpoint_finite"] != "1":
                        raise ValueError(f"Latency call validation failed: {path}")
                    if r["fixture_id"] != f"noise{int(r['block_id']) % 20:02d}":
                        raise ValueError(f"Frozen noise cycle mismatch: {path}")
                    wall = finite_number(float(r["wall_ms"]), "complete reconstruction wall", True)
                    start, ready, checked = (int(r[k]) for k in
                                             ("t_start_ns", "t_result_ready_ns", "t_light_check_end_ns"))
                    if not start < ready <= checked:
                        raise ValueError(f"Nonmonotonic latency call endpoints: {path}")
                    close(wall, (ready - start) / 1e6, "wall and CLOCK_MONOTONIC endpoints")
            x = [float(r["wall_ms"]) for r in measured]
            values[method].extend(x)
            sessions[method].append({"session": session, "N": 200, "mean_ms": statistics.mean(x)})
            sources.append(path)
    if any(len(x) != 600 for x in values.values()) or sum(map(len, values.values())) != 1800:
        raise ValueError("Complete 600 measured calls per method required")
    return {"values": values, "sessions": sessions, "sources": sources,
            "excluded_protocol_warmups": 450, "measured_calls": 1800}


def type7(values, q):
    x = sorted(values)
    h = (len(x) - 1) * q
    lower, upper = math.floor(h), math.ceil(h)
    return x[lower] + (h - lower) * (x[upper] - x[lower])


def latency_statistics(latency):
    return {
        m: {
            "N": len(x), "mean_ms": statistics.mean(x),
            "sample_SD_ms": statistics.stdev(x), "median_ms": statistics.median(x),
            "P95_type7_ms": type7(x, .95), "min_ms": min(x), "max_ms": max(x),
            "sessions": latency["sessions"][m],
        }
        for m, x in latency["values"].items()
    }


def style():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
        "axes.linewidth": .8, "lines.linewidth": 1.0,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.unicode_minus": True, "figure.dpi": 120, "savefig.dpi": 600,
        "svg.fonttype": "none", "svg.hashsalt": "sonar-formal-power-20261010",
        "path.simplify": False, "agg.path.chunksize": 0,
    })


def export(fig, out, stem, description):
    png = out / f"{stem}.png"
    svg = out / f"{stem}.svg"
    fig.savefig(png, format="png", dpi=600,
                metadata={"Title": stem, "Description": description, "Software": "plot_power_results.py"})
    fig.savefig(svg, format="svg",
                metadata={"Title": stem, "Description": description, "Creator": "plot_power_results.py", "Date": None})
    return [png, svg]


def summary_plot(windows, out):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.2))
    fig.subplots_adjust(left=.105, right=.99, bottom=.23, top=.82, wspace=.37)
    summaries = {}
    fields = ("average_power_W", "energy_per_block_J")
    ylabels = ("Average module power (W)", "Module energy per block (J/block)")
    for m in METHODS:
        summaries[m] = {}
        for field in fields:
            x = [w[field] for w in windows if w["method"] == m]
            if len(x) != 3:
                raise ValueError("Exactly three independent windows per method are required")
            summaries[m][field] = {"values": x, "mean": statistics.mean(x), "sample_SD": statistics.stdev(x)}
    for panel, (ax, field, ylabel) in enumerate(zip(axes, fields, ylabels)):
        maximum = 0.
        for i, method in enumerate(METHODS):
            s = summaries[method][field]
            if s["mean"] - s["sample_SD"] < 0:
                raise ValueError("Mean-SD crosses zero; redesign the interval display before plotting")
            positions = i + np.array([-.11, 0., .11])
            ax.scatter(positions, s["values"], s=34, marker=MARKERS[method],
                       facecolor="white", edgecolor=COLORS[method], linewidth=1.1, zorder=4)
            ax.errorbar(i, s["mean"], yerr=s["sample_SD"], fmt="_", color="#222222",
                        markersize=13, capsize=4, capthick=1., linewidth=1., zorder=5)
            maximum = max(maximum, max(s["values"]), s["mean"] + s["sample_SD"])
        ax.set_xticks(range(3), [LABELS[m] for m in METHODS])
        ax.set_xlim(-.45, 2.45)
        ax.set_ylim(0, maximum * 1.18)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="#D9D9D9", linewidth=.5)
        ax.set_axisbelow(True)
        ax.text(.02, .96, f"({chr(97 + panel)})", transform=ax.transAxes,
                va="top", ha="left", fontsize=11)
    legend = [
        Line2D([], [], marker="o", color="#555555", markerfacecolor="white", linestyle="none", markersize=5, label="Individual window (n = 3)"),
        Line2D([], [], marker="_", color="#222222", markersize=11, linestyle="-", linewidth=1., label="Mean ± sample SD"),
    ]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.54, .99), ncol=2, frameon=False)
    outputs = export(fig, out, "module_power_energy",
                     "Three independent windows per method; mean and sample SD. Module VDD_IN, no idle subtraction.")
    plt.close(fig)
    return outputs, summaries


def jetson_summary_plot(latency, summaries, out):
    import matplotlib.pyplot as plt
    import numpy as np
    fig = plt.figure(figsize=(6.8, 5.0))
    grid = fig.add_gridspec(2, 2, left=.105, right=.985, bottom=.18, top=.91,
                           hspace=.7, wspace=.37, height_ratios=(1.05, 1.0))
    top = fig.add_subplot(grid[0, :])
    bottom = [fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    data = [latency["values"][m] for m in METHODS]
    # Deterministic offsets make every measured call a visible input to the
    # scatter; no outlier trimming, winsorization, or thinning is applied.
    for i, method in enumerate(METHODS):
        x = np.array(data[i])
        jitter = ((np.arange(len(x)) * 37) % 101 - 50) / 50 * .15
        top.scatter(i + jitter, x, s=5, color=COLORS[method], marker=MARKERS[method],
                    alpha=.17, linewidth=0, zorder=1, rasterized=False)
    boxes = top.boxplot(data, positions=range(3), widths=.42, whis=(0, 100), showfliers=True,
                        patch_artist=True, showmeans=True,
                        medianprops={"color": "#111111", "linewidth": 1.1},
                        meanprops={"marker": "_", "markeredgecolor": "#111111",
                                   "markerfacecolor": "#111111", "markersize": 10},
                        whiskerprops={"color": "#444444", "linewidth": .8},
                        capprops={"color": "#444444", "linewidth": .8})
    for box, method in zip(boxes["boxes"], METHODS):
        box.set(facecolor=COLORS[method] + "35", edgecolor=COLORS[method], linewidth=1.)
    top.set_xticks(range(3), [LABELS[m] for m in METHODS])
    top.set_xlim(-.55, 2.55)
    top.set_ylim(0, max(map(max, data)) * 1.13)
    top.set_ylabel("Full reconstruction time (ms)")
    top.grid(axis="y", color="#D9D9D9", linewidth=.5)
    top.set_axisbelow(True)
    top.text(.02, .96, "(a)", transform=top.transAxes, va="top", ha="left", fontsize=11)
    top.text(.99, 1.04, "600 calls/method; whiskers: min–max",
             transform=top.transAxes, ha="right", fontsize=9)
    for panel, (ax, field, ylabel) in enumerate(zip(
            bottom, ("average_power_W", "energy_per_block_J"),
            ("Average module power (W)", "Module energy per block (J/block)")), start=1):
        maximum = 0.
        for i, method in enumerate(METHODS):
            s = summaries[method][field]
            ax.scatter(i + np.array([-.11, 0., .11]), s["values"], s=32, marker=MARKERS[method],
                       facecolor="white", edgecolor=COLORS[method], linewidth=1.1, zorder=4)
            ax.errorbar(i, s["mean"], yerr=s["sample_SD"], fmt="_", color="#222222",
                        markersize=12, capsize=4, capthick=1., linewidth=1., zorder=5)
            maximum = max(maximum, max(s["values"]), s["mean"] + s["sample_SD"])
        ax.set_xticks(range(3), [LABELS[m] for m in METHODS])
        ax.set_xlim(-.45, 2.45)
        ax.set_ylim(0, maximum * 1.18)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="#D9D9D9", linewidth=.5)
        ax.set_axisbelow(True)
        ax.text(.02, .96, f"({chr(97 + panel)})", transform=ax.transAxes,
                va="top", ha="left", fontsize=11)
    fig.text(.5, .09, "(b,c): n = 3 windows/method; mean ± sample SD", ha="center", fontsize=9)
    outputs = export(fig, out, "jetson_summary",
                     "Same campaign sources: all 600 complete measured calls per method; 3 independent module-power windows per method.")
    plt.close(fig)
    return outputs


def trace_plot(windows, out, temperature=False):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator
    by_key = {(w["session"], w["method"]): w for w in windows}
    field = "max_temperature_C" if temperature else "power_W"
    series = "temperature_rows" if temperature else "power_rows"
    if temperature and any(not w[series] for w in windows):
        raise ValueError("Temperature graph requires available readouts in all nine windows")
    all_times = [(r["timestamp_ns"] - w["t0_ns"]) / 1e9 for w in windows for r in w[series]]
    all_values = [r[field] for w in windows for r in w[series]]
    x_low, x_high = math.floor(min(all_times) / 10) * 10, math.ceil(max(all_times) / 10) * 10
    y_high = max(85., max(all_values) * 1.06) if temperature else max(all_values) * 1.08
    fig, axes = plt.subplots(3, 3, figsize=(6.8, 6.1), sharex=True, sharey=True)
    fig.subplots_adjust(left=.13, right=.985, bottom=.105, top=.88, hspace=.19, wspace=.14)
    for row, session in enumerate(SESSIONS):
        for column, method in enumerate(METHODS):
            ax = axes[row, column]
            w = by_key[session, method]
            idle = (w["idle_ns"] - w["t0_ns"]) / 1e9
            warm = (w["warm_ns"] - w["t0_ns"]) / 1e9
            end = (w["t1_ns"] - w["t0_ns"]) / 1e9
            ax.axvspan(idle, warm, facecolor=PHASE_COLORS["Idle"], edgecolor="none", zorder=0)
            ax.axvspan(warm, 0., facecolor=PHASE_COLORS["Warm-up"], edgecolor="none", zorder=0)
            ax.axvspan(0., end, facecolor=PHASE_COLORS["Measured window"], edgecolor="none", zorder=0)
            x = np.array([(r["timestamp_ns"] - w["t0_ns"]) / 1e9 for r in w[series]])
            y = np.array([r[field] for r in w[series]])
            # Every valid original row is a line vertex; markers only distinguish
            # methods in grayscale and do not resample or smooth the data.
            ax.plot(x, y, color=COLORS[method], linestyle=LINESTYLES[method], linewidth=.8,
                    marker=MARKERS[method], markersize=2.5, markerfacecolor="white",
                    markeredgewidth=.6, markevery=max(1, len(x) // 24), zorder=2)
            ax.axvline(0., color="#444444", linewidth=.65, linestyle="-", zorder=1)
            if temperature:
                ax.axhline(80., color="#555555", linewidth=.8, linestyle="--", zorder=3)
            ax.set_xlim(x_low, x_high)
            ax.set_ylim(0., y_high)
            ax.set_xticks([-90, -60, 0, 60, 120])
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
            ax.grid(axis="y", color="#D9D9D9", linewidth=.4)
            ax.set_axisbelow(True)
            if row == 0:
                ax.set_title(LABELS[method], pad=6)
        position = axes[row, 0].get_position()
        fig.text(.065, (position.y0 + position.y1) / 2,
                 f"Session {session + 1}", va="center", ha="center", rotation=90, fontsize=9)
    fig.supxlabel("Time relative to measurement start (s)", fontsize=10, y=.025)
    fig.supylabel("Maximum available sensor temperature (°C)" if temperature else "Module VDD_IN power (W)",
                  fontsize=10, x=.007)
    handles = [Patch(facecolor=color, edgecolor="#999999", linewidth=.4, label=name)
               for name, color in PHASE_COLORS.items()]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.55, .985), ncol=3, frameon=False)
    stem = "temperature_time_all_windows" if temperature else "power_time_all_windows"
    outputs = export(fig, out, stem, "Complete valid original telemetry reads, no smoothing; shared axes and actual runner phase boundaries.")
    plt.close(fig)
    return outputs


def captions(windows, summaries, include_temperature, latency_stats):
    power_means = [summaries[m]["average_power_W"]["mean"] for m in METHODS]
    energy_means = [summaries[m]["energy_per_block_J"]["mean"] for m in METHODS]
    result = {
        "jetson_summary": {
            "english": (
                "The three-panel summary reports complete reconstruction latency and module power/energy from the same frozen three-method campaign. "
                "Panel (a) retains all 600 measured calls per method across three sessions; the box marks median and interquartile range, the central short mark shows the mean, and whiskers span the actual minimum and maximum with no outlier removal. "
                "Only the 50 predefined warm-up calls per session are excluded from the measured distribution, and the underlying wall-time CSV files also supply the corresponding table statistics. "
                "Panels (b,c) show three measurement windows per method, with individual-window symbols and mean ± sample standard deviation; energy uses the actual VDD_IN submission-and-drain interval without idle subtraction. "
                "Dense LFM/CDM-ZC use cached full-rate processing and FDM-LFM uses D4/sinc16 on the same output grid; the timings describe reconstruction, not acoustic frame rates."
            ),
            "chinese": (
                "三面板图从同一冻结的三方法测试数据中呈现完整重建时延及模块功率/能耗。"
                "(a)保留每种方法三个会话中的全部600次正式调用：箱体为中位数与四分位范围，短标记为均值，须线覆盖实际最小与最大值，不删除离群样本。"
                "仅按预定定义排除每会话50次暖机调用；对应表格统计同样来自这些原始wall时间CSV。"
                "(b,c)每种方法为三个独立测量窗口的散点、均值及样本标准差，能耗按实际VDD_IN提交及末块完成窗口积分，不扣空闲功率。"
                "Dense LFM/CDM-ZC采用缓存全速处理，FDM-LFM采用D4/sinc16并输出相同网格；时延指重建计算，不是声学帧率。"
            ),
            "latency_statistics_from_same_sources": latency_stats,
        },
        "module_power_energy": {
            "english": (
                f"Across three independent windows per method, module-average power spans {min(power_means):.2f}–{max(power_means):.2f} W and energy per reconstructed block spans {min(energy_means):.4f}–{max(energy_means):.4f} J. "
                "Symbols show individual windows (n = 3 per method), and black marks show the mean ± sample standard deviation across windows, not across telemetry samples. "
                "Power is the current Jetson module VDD_IN reading, and energy is integrated over each actual submission-and-drain window without idle subtraction. "
                "Dense LFM and CDM-ZC use cached full-rate reconstruction; FDM-LFM uses D4 decimation with 16-tap local sinc interpolation on the same 321 × 401 output grid. "
                "These are module-level reconstruction measurements, not carrier-board or complete-sonar energy measurements."
            ),
            "chinese": (
                f"每种方法的三个独立窗口中，平均模块功率为{min(power_means):.2f}–{max(power_means):.2f} W，每个重建块的模块能耗为{min(energy_means):.4f}–{max(energy_means):.4f} J。"
                "散点为各独立窗口（每种方法n=3），黑色标记为窗口间均值及样本标准差，而不是遥测采样点间统计。"
                "使用Jetson模块VDD_IN当前功率，按实际提交及末块完成窗口积分，不扣除空闲功率。"
                "Dense LFM和CDM-ZC采用缓存全速处理，FDM-LFM采用D4抽取与16-tap局部sinc插值，共同输出321×401网格。"
                "统计范围为模块重建计算，不代表载板或完整声呐能耗。"
            ),
        },
        "power_time_all_windows": {
            "english": (
                "The complete telemetry traces retain the idle, warm-up and measured-window transitions for all nine runs. "
                "Columns identify the reconstruction method and rows identify the session; time zero is the actual measured-window start, and shaded intervals follow the runner's monotonic phase timestamps. "
                "All valid current VDD_IN readings are retained without smoothing or downsampling, and both axes are shared across panels. "
                "Repeated readings are permitted and do not imply independent 100-ms sensor measurements."
            ),
            "chinese": (
                "九次运行的完整遥测曲线保留了空闲、预热及测量窗口的阶段变化。"
                "列对应处理方法，行对应独立会话；零时刻为实际测量窗口起点，背景区域按runner单调时钟阶段事件划分。"
                "所有有效VDD_IN当前读值均保留，不平滑、不抽样，各面板共享坐标范围。"
                "允许重复读值，100 ms采集请求不等同于独立传感器测量带宽。"
            ),
        },
    }
    if include_temperature:
        t = [r["max_temperature_C"] for w in windows for r in w["temperature_rows"]]
        result["temperature_time_all_windows"] = {
            "english": (
                f"The maximum among available sensor readouts ranges from {min(t):.1f} to {max(t):.1f} °C over the logged phases. "
                "All available readouts are retained on the same session/method layout and phase boundaries as the power traces; the dashed line marks the 80 °C safety-stop threshold. "
                "This available-sensor maximum is not a dedicated GPU-temperature measurement and does not establish constant clock frequency or the absence of throttling."
            ),
            "chinese": (
                f"记录阶段中，可获得传感器的最高温度读值为{min(t):.1f}–{max(t):.1f} °C。"
                "全部可用读值按照与功率图相同的会话、方法及阶段边界绘制，虚线为80 °C安全停止阈值。"
                "该量为可获得传感器的最高读值，不专指GPU温度，也不能据此推断恒定频率或无节流。"
            ),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, help="Default: repository outputs/jetson/diagnostics; refuse existing outputs.")
    parser.add_argument("--temperature", action="store_true", help="Also require and plot temperature in all nine windows.")
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    windows = load_campaign(run_root)  # No placeholders or partial-campaign plots.
    latency = load_latency_campaign(run_root)
    latency_stats = latency_statistics(latency)
    stems = list(OUTPUT_STEMS) + (["temperature_time_all_windows"] if args.temperature else [])
    out = args.out.resolve() if args.out else RELEASE_ROOT / "outputs" / "jetson" / "diagnostics"
    expected = [out / f"{stem}.{suffix}" for stem in stems for suffix in ("png", "svg")]
    expected += [out / "captions.json", out / "power_window_metrics.csv", out / "jetson_summary_latency_stats.csv", out / "FIGURE_MANIFEST.json"]
    if any(p.exists() for p in expected):
        raise FileExistsError("Refusing to overwrite generated figure artifacts; choose a fresh --out directory")
    if args.temperature and any(not w["temperature_rows"] for w in windows):
        raise ValueError("Temperature requested but not available in all nine windows")
    style()
    out.mkdir(parents=True, exist_ok=True)
    outputs, summaries = summary_plot(windows, out)
    outputs += jetson_summary_plot(latency, summaries, out)
    outputs += trace_plot(windows, out)
    if args.temperature:
        outputs += trace_plot(windows, out, temperature=True)
    caption_path = out / "captions.json"
    caption_path.write_text(json.dumps(captions(windows, summaries, args.temperature, latency_stats), indent=2,
                                       ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    metrics_path = out / "power_window_metrics.csv"
    with metrics_path.open("x", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["session", "method", "average_power_W", "energy_per_block_J",
                                              "energy_J", "duration_s", "completed_blocks"])
        writer.writeheader()
        for w in windows:
            writer.writerow({k: w[k] for k in writer.fieldnames})
    latency_path = out / "jetson_summary_latency_stats.csv"
    with latency_path.open("x", newline="", encoding="utf-8") as f:
        fields = ["method", "N", "mean_ms", "sample_SD_ms", "median_ms", "P95_type7_ms", "min_ms", "max_ms"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for method in METHODS:
            writer.writerow({"method": method, **{k: latency_stats[method][k] for k in fields if k != "method"}})
    import matplotlib
    import numpy
    from matplotlib.font_manager import FontProperties, findfont
    font_path = Path(findfont(FontProperties(family=["Times New Roman", "DejaVu Serif"])))
    sources = sorted({p for w in windows for p in w["sources"]} | set(latency["sources"]))
    manifest = {
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "run_id": run_root.name, "independent_windows_per_method": 3, "window_count": 9,
        "statistics": "mean and sample SD (ddof=1) across independent windows",
        "method_labels": LABELS, "FDM_engineering_chain": "D4, fs_low=100 kHz, local sinc16",
        "module_power_scope": "Current VDD_IN, not carrier/whole-sonar; primary energy has no idle subtraction",
        "raw_read_policy": "Every valid row is retained in original order; path.simplify=False; no smoothing or downsampling",
        "marker_policy": "Sparse decorative method markers do not alter retained line vertices",
        "dimensions_inches": {"jetson_summary": [6.8, 5.0], "module_power_energy": [6.8, 3.2], "traces": [6.8, 6.1]},
        "font_sizes_pt": {"labels": 10, "ticks": 9, "legend": 9, "panel_letters": 11},
        "export": {"png_dpi": 600, "svg_text": True, "embedded_raster": False, "PDF": False},
        "runtime": {"python": sys.version, "matplotlib": matplotlib.__version__, "numpy": numpy.__version__, "platform": platform.platform()},
        "resolved_font": {"file": font_path.name, "sha256": sha(font_path)},
        "auxiliary_reader_errors_by_window": [
            {"session": w["session"], "method": w["method"], "errors": w["auxiliary_reader_errors"]}
            for w in windows if w["auxiliary_reader_errors"]
        ],
        "sources": [{"path_relative_to_run_root": p.relative_to(run_root).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)} for p in sources],
        "outputs": [{"path_relative_to_output_root": p.relative_to(out).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)}
                    for p in outputs + [caption_path, metrics_path, latency_path]],
        "window_summary": summaries,
        "latency_statistics": latency_stats,
        "latency_data_policy": {"formal_calls": 1800, "calls_per_method": 600, "protocol_warmups_excluded": 450,
                                "outlier_removal": False, "whiskers": "actual min/max", "scatter": "all measured calls"},
        "post_export_visual_QA": "Required from primary agent; not asserted by generator.",
    }
    manifest_path = out / "FIGURE_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"figures_ready": True, "out": str(out), "windows": 9,
                      "figure_count": len(stems), "visual_QA_pending": True}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
