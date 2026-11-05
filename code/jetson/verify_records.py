"""Recalculate the paper's Jetson table from released records, without a GPU."""
from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

CODE = Path(__file__).resolve().parent
RELEASE = CODE.parents[1]
DEFAULT_RECORDS = RELEASE / "data" / "jetson"
METHODS = ("dense", "cdm", "fdm")
POWER = re.compile(r"\bVDD_IN\s+(\d+)mW/(\d+)mW")
TEMP = re.compile(r"[A-Za-z0-9_]+@([0-9.]+)C")
FROZEN_SOURCES = (
    "formal_collection_20261010/core/CMakeLists.txt",
    "formal_collection_20261010/core/main.cpp",
    "formal_collection_20261010/telemetry.py",
    "implementation_v2/core/engine.cu",
    "implementation_v2/core/engine.h",
    "implementation_v3/core/engine_v3.cu",
    "implementation_v3/core/engine_v3.h",
)
MODEL_SHA256 = "4b280863495d4d0f677f36dac79736cedd9d2ba009a3cacf41849cfce7db8f64"


def require(test, message):
    if not test:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def close(a, b, label, tolerance=1e-11):
    require(math.isclose(float(a), float(b), abs_tol=tolerance, rel_tol=tolerance),
            f"Inconsistent {label}: {a} vs {b}")


def type7(values, q):
    x = sorted(values)
    h = (len(x) - 1) * q
    lo, hi = math.floor(h), math.ceil(h)
    return x[lo] + (h - lo) * (x[hi] - x[lo])


def latency_values(path, method, session):
    rows = read_csv(path)
    require(len(rows) == 250, f"Wrong total row count: {path.name}")
    result = []
    for phase, count in (("warmup", 50), ("measure", 200)):
        group = [r for r in rows if r["phase"] == phase]
        require([int(r["block_id"]) for r in group] == list(range(count)),
                f"Incomplete {method}/{session}/{phase}")
        for r in group:
            require(r["method"] == method and int(r["session"]) == session and
                    r["status"] == "ok" and r["shape_pass"] == "1" and r["checkpoint_finite"] == "1",
                    f"Failed or mislabeled latency call: {method}/{session}/{phase}")
            require(r["fixture_id"] == f"noise{int(r['block_id']) % 20:02d}", "Wrong noise cycle")
            start, ready, checked = (int(r[k]) for k in
                                     ("t_start_ns", "t_result_ready_ns", "t_light_check_end_ns"))
            require(start < ready <= checked, "Nonmonotonic latency endpoints")
            wall = float(r["wall_ms"])
            close(wall, (ready-start)/1e6, "wall and monotonic endpoints")
            if phase == "measure":
                result.append(wall)
    return result


def power_metrics(directory, method, session):
    window = read_json(directory/"runner"/"window.json")
    logged = read_json(directory/"POWER_INTEGRAL.json")
    invocation = read_json(directory/"invocation.json")
    require(window["method"] == method and window["session"] == session and
            invocation["returncode"] == 0 and not invocation["whole_runner_timeout"] and
            window["workload_valid"] is True and logged["valid"] is True and
            logged["idle_subtracted"] is False, f"Invalid power window: {method}/{session}")
    t0, t1, n = (window[k] for k in ("t0_ns", "t1_ns", "completed_blocks"))
    require(n > 0 and t1-t0 >= 120_000_000_000 and all(window[k] == logged[k] for k in
            ("t0_ns", "t1_ns", "completed_blocks")), "Power integral/window boundary mismatch")
    phase_counts = {}
    with (directory/"runner"/"power_blocks.csv").open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            require(r["status"] == "ok" and r["method"] == method and int(r["session"]) == session,
                    "Invalid power-campaign block")
            phase_counts[r["phase"]] = phase_counts.get(r["phase"], 0) + 1
            if r["phase"] == "power_measure":
                require(t0 <= int(r["t_start_ns"]) < int(r["t_result_ready_ns"]) <=
                        int(r["t_light_check_end_ns"]) <= t1, "Power block lies outside window")
    require(set(phase_counts) == {"power_warmup", "power_measure"} and
            phase_counts.get("power_measure") == n and
            phase_counts.get("power_warmup") == window["warmup_completed_blocks"], "Wrong completed block count")
    events = {e["event"]: e for e in invocation["events"]}
    require(events["idle_begin"]["monotonic_ns"] < events["warmup_begin"]["monotonic_ns"] < t0 < t1 and
            events["window_start"]["t0_ns"] == t0 and events["window_end"]["t1_ns"] == t1,
            "Inconsistent power phase event")
    with (directory/"power_raw.jsonl").open(encoding="utf-8") as f:
        raw = [json.loads(line) for line in f]
    values = []
    for row in raw:
        match = POWER.search(row["line"])
        temps = [float(t) for t in TEMP.findall(row["line"])]
        values.append((row["timestamp_ns"], int(match[1])/1000 if match else None,
                       max(temps) if temps else None))
    require(all(b[0] > a[0] for a,b in zip(values, values[1:])), "Non-increasing telemetry timestamps")
    before = [i for i,r in enumerate(values) if r[0] <= t0 and r[1] is not None]
    after = [i for i,r in enumerate(values) if r[0] >= t1 and r[1] is not None]
    require(before and after, "Unbracketed telemetry boundaries")
    region = values[before[-1]:after[0]+1]
    require(all(r[1] is not None and (r[2] is None or r[2] < 80) for r in region),
            "Missing power or safety threshold reached")
    gap_ns = max(b[0]-a[0] for a,b in zip(region,region[1:]))
    require(gap_ns <= 500_000_000, "Power telemetry gap above 500 ms")
    times = [r[0] for r in region]

    def boundary(t):
        j = bisect.bisect_left(times, t)
        if times[j] == t:
            return t, region[j][1]
        (a,x,_),(b,y,_) = region[j-1:j+1]
        return t, x+(t-a)/(b-a)*(y-x)

    points = [boundary(t0)] + [(r[0],r[1]) for r in region if t0<r[0]<t1] + [boundary(t1)]
    energy = math.fsum((b[0]-a[0])*(a[1]+b[1])*.5e-9 for a,b in zip(points,points[1:]))
    duration = (t1-t0)*1e-9
    average = energy/duration
    per_block = energy/n
    close(energy, logged["energy_J"], "integrated energy")
    close(average, logged["average_module_power_W"], "window average power")
    close(per_block, logged["energy_J_per_block"], "window energy per block")
    return {"power_W":average,"energy_J_per_block":per_block,"energy_J":energy,
            "completed_blocks":n,"duration_s":duration,"maximum_interval_ms":gap_ns/1e6,
            "maximum_bracket_temperature_C":max((r[2] for r in region if r[2] is not None),default=None)}


def memory_metrics(directory):
    m = read_json(directory/"runner"/"memory.json")
    ledger = sum(int(r["bytes"]) for r in read_csv(directory/"runner"/"buffers.csv"))
    require(ledger == m["explicit_peak_bytes"] == m["explicit_live_bytes"] and
            m["explicit_peak_bytes"] == m["non_workspace_bytes"]+m["shared_plan_workspace_bytes"],
            "Inconsistent explicit CUDA allocation ledger")
    require(m["not_total_module_peak"] and m["unified_memory_quantities_must_not_be_summed"],
            "Incorrect memory scope")
    samples = read_json(directory/"power_samples.json")["rows"]
    ram = [r["system_RAM_used_bytes"] for r in samples if r.get("system_RAM_used_bytes") is not None]
    ram_report = read_json(directory/"SYSTEM_RAM_SAMPLED.json")
    require(len(ram)==ram_report["sample_count"] and max(ram,default=None)==ram_report["sampled_max_used_bytes"],
            "Sampled system RAM receipt mismatch")
    return m["explicit_peak_bytes"] / 2**20


def verify(records=DEFAULT_RECORDS):
    manifest_path = records/"provenance"/"FROZEN_FORMAL_FINAL.json"
    freeze = read_json(manifest_path)
    frozen = {r["path"]:r for r in freeze["files"]}
    for name in FROZEN_SOURCES:
        require((CODE/name).stat().st_size == frozen[name]["bytes"] and
                sha(CODE/name) == frozen[name]["sha256"], f"Frozen source mismatch: {name}")
    require(sha(CODE/"reference_snapshot"/"numerical_model.py") == MODEL_SHA256,
            "Frozen fixture model source hash changed")
    validation = read_json(records/"provenance"/"IMPLEMENTATION_RECEIPT_FORMAL.json")
    require(validation["binary_sha256"] == freeze["binary_sha256"] and
            validation["frozen_manifest_sha256"] == sha(manifest_path) and
            validation["all_passed"] and validation["comparison_count"] == 75 and
            validation["equivalence_claimed"] is False, "Historical 75-case validation receipt mismatch")
    latency_summary = read_json(records/"summary"/"SUMMARY_FORMAL.json")["latency"]
    power_summary = read_json(records/"summary"/"POWER_SUMMARY.json")["methods"]
    figure_latency = {r["method"]:r for r in read_csv(records/"summary"/"jetson_summary_latency_stats.csv")}
    figure_power = {(r["method"],int(r["session"])):r for r in read_csv(records/"summary"/"power_window_metrics.csv")}
    require(set(figure_power) == {(m,s) for m in METHODS for s in range(3)}, "Incomplete figure power summary")
    methods = {}
    for method in METHODS:
        times, windows, memory = [], [], []
        for session in range(3):
            key = Path(f"session{session}")/method
            times.extend(latency_values(records/"latency"/key/"runner"/"latency_raw.csv",method,session))
            w = power_metrics(records/"power"/key,method,session)
            windows.append(w)
            f = figure_power[(method,session)]
            close(w["power_W"],f["average_power_W"],"figure power")
            close(w["energy_J_per_block"],f["energy_per_block_J"],"figure energy")
            require(w["completed_blocks"]==int(f["completed_blocks"]),"Figure block count mismatch")
            memory.append(memory_metrics(records/"memory"/key))
        require(len(times)==600 and len(set(memory))==1,"Incomplete timing or changing explicit CUDA allocation")
        p = [w["power_W"] for w in windows]
        e = [w["energy_J_per_block"] for w in windows]
        row = {"N":len(times),"mean_ms":statistics.mean(times),"sample_SD_ms":statistics.stdev(times),
               "median_ms":statistics.median(times),"P95_type7_ms":type7(times,.95),
               "min_ms":min(times),"max_ms":max(times),"explicit_CUDA_MiB":memory[0],
               "power_windows":3,"power_W_mean":statistics.mean(p),"power_W_sample_SD":statistics.stdev(p),
               "energy_J_per_block_mean":statistics.mean(e),"energy_J_per_block_sample_SD":statistics.stdev(e),
               "completed_power_blocks":sum(w["completed_blocks"] for w in windows),"windows":windows}
        require(int(figure_latency[method]["N"])==len(times), "Wrong figure call count")
        for k in ("mean_ms","sample_SD_ms","median_ms","P95_type7_ms","min_ms","max_ms"):
            close(row[k],figure_latency[method][k],f"figure {method}/{k}")
        for k,source_key in (("mean_ms","mean_ms"),("sample_SD_ms","sample_SD_ms"),
                             ("median_ms","median_ms"),("P95_type7_ms","P95_ms"),("max_ms","max_ms")):
            close(row[k],latency_summary[method][source_key],f"latency summary {method}/{k}")
        for prefix,summary_key in (("power_W","average_module_power_W_window_mean_and_sample_SD"),
                                   ("energy_J_per_block","energy_J_per_block_window_mean_and_sample_SD")):
            for stat in ("mean","sample_SD"):
                close(row[f"{prefix}_{stat}"],power_summary[method][summary_key][stat],f"power summary {method}/{prefix}/{stat}")
        methods[method] = row
    return {"passed":True,"mode":"offline reanalysis; no new device measurements",
            "frozen_sources_verified":len(FROZEN_SOURCES),"historical_own_reference_checks":75,
            "latency_calls":sum(v["N"] for v in methods.values()),"power_windows":9,
            "completed_power_blocks":sum(v["completed_power_blocks"] for v in methods.values()),
            "historical_binary_sha256":freeze["binary_sha256"],"methods":methods}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records",type=Path,default=DEFAULT_RECORDS)
    args=parser.parse_args()
    print(json.dumps(verify(args.records),indent=2,allow_nan=False))


if __name__=="__main__":
    main()
