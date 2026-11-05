"""Verify the release manifest and selected data-to-result contracts."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
MANIFEST = ROOT / "SHA256SUMS.txt"
GEOMETRY = DATA / "shipwreck" / "geometry"
GEOMETRY_MANIFEST = GEOMETRY.parent / "GEOMETRY_SHA256SUMS.txt"
DESKTOP = DATA / "desktop-timing" / "overlay"
DESKTOP_MANIFEST = DESKTOP.parent / "HDF5_SHA256SUMS.txt"
METHODS = {"dense", "cdm", "fdm"}
IGNORED_PARTS = {".git", ".venv", "outputs", "build", "__pycache__", ".pytest_cache"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def files():
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if (path.is_file() and path != MANIFEST
                and not any(path.is_relative_to(folder) for folder in (GEOMETRY, DESKTOP))
                and not IGNORED_PARTS.intersection(relative.parts)):
            yield path


def verify_asset(directory, manifest, archive):
    require(directory.is_dir(), f"Extract {archive} at the repository root first.")
    expected = dict((name, hash_) for hash_, name in (
        line.split("  ", 1) for line in manifest.read_text(encoding="utf-8").splitlines()))
    actual = {p.relative_to(directory).as_posix(): digest(p)
              for p in directory.rglob("*") if p.is_file()}
    require(actual == expected, f"{archive}: checksum or file set mismatch.")
    return len(actual)


def rows(name):
    with (DATA / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def close(a, b, message):
    require(np.isclose(a, b, rtol=1e-12, atol=1e-12), message)


def trial_groups(folder, metrics, expected_groups):
    trials = rows(f"{folder}/all_trials.csv")
    summaries = rows(f"{folder}/summary_mean_ci.csv")
    index = ("noise_normalization", "nominal_snr_db", "method")
    groups, pairs = defaultdict(list), defaultdict(list)
    for row in trials:
        key = tuple(row[k] for k in index)
        groups[key].append(row)
        pairs[(*key[:2], row.get("phase_repetition", ""), row["noise_repetition"])].append(row)
    require(len(groups) == len(summaries) == expected_groups, f"{folder}: summary group count")
    for pair in pairs.values():
        require(len(pair) == 3 and {q["method"] for q in pair} == METHODS, f"{folder}: pairing")
        require(len({q["noise_seed"] for q in pair}) == 1, f"{folder}: paired seed")
    for s in summaries:
        g = groups[tuple(s[k] for k in index)]
        require(len(g) == int(s["replications"]) == 30, f"{folder}: trial count")
        for key in metrics:
            x = np.array([float(q[key]) for q in g])
            close(float(s[key]), x.mean(), f"{folder}: {key} mean")
            close(float(s[f"{key}_ci95_halfwidth"]), t.ppf(.975, 29) * x.std(ddof=1) / np.sqrt(30),
                  f"{folder}: {key} interval")
    return len(trials)


def science():
    report = {}
    report["point_trials"] = trial_groups("point30m", [
        "cnr_linear_magnitude", "noise_induced_magnitude_NMSE_db", "detection_probability",
        "false_detections"], 24)
    report["dense_trials"] = trial_groups("dense1000x50", [
        "noise_induced_magnitude_NMSE_db", "dense_magnitude_correlation",
        "truth_magnitude_correlation"], 24)
    report["shipwreck_trials"] = trial_groups("shipwreck", [
        "dense_magnitude_correlation", "range_profile_correlation_to_clean_Dense",
        "angle_profile_correlation_to_clean_Dense", "support_leakage"], 6)

    with np.load(DATA / "point30m/point30m_main_minus15dB_arrays.npz", allow_pickle=False) as d:
        r, a = d["ranges_m"], d["angles_deg"]
        require(len(d["true_ranges_m"]) == len(d["true_angles_deg"]) == 18, "Point target count")
        require(30 < d["true_ranges_m"].min() <= d["true_ranges_m"].max() < 31, "Point range")
        for method in METHODS:
            x = d[method]
            require(x.shape == (r.size, a.size) and np.isfinite(x).all() and np.iscomplexobj(x),
                    f"Figure 10 shape/type: {method}")
        report["figure10_shape"] = [r.size, a.size]

    density = [r for r in rows("density/metrics.csv") if r["factor"] in {"scatter_count", "range_layers"}]
    groups = defaultdict(list)
    for row in density:
        groups[(row["factor"], row["factor_value"])].append(row)
    require(len(density) == 99 and len(groups) == 11, "Figure 12 selection")
    for group in groups.values():
        require(len({r["scene_seed"] for r in group}) == 3 and len(group) == 9, "Figure 12 repeat count")
        require(all(float(r["background_fraction"]) >= .7202 for r in group), "Figure 12 guard")
    report["figure12_rows"] = len(density)

    c = next(r for r in rows("coverage/dense_final_audit.csv") if r["design"] == "implemented")
    require(int(c["sector_count"]) == 6, "Implemented sector count")
    require(-6.6 <= float(c["dense_worst_coverage_db"]) < -6.5, "Coverage condition")
    require(-24.4 < float(c["dense_worst_in_fov_grating_db"]) <= -24, "Grating condition")
    report["coverage_and_grating_db"] = [float(c["dense_worst_coverage_db"]),
                                        float(c["dense_worst_in_fov_grating_db"])]

    controlled = rows("das_omp/controlled_metrics.csv")
    groups = defaultdict(list)
    for row in controlled:
        groups[(row["range_step_m"], row["angle_step_deg"], row["noise_replication"])].append(row)
    require(len(controlled) == 12 and len(groups) == 6, "DAS/OMP comparison count")
    for pair in groups.values():
        require({r["method"] for r in pair} == {"DAS", "OMP"}
                and len({r["measurement_sha256"] for r in pair}) == 1, "DAS/OMP matched measurements")
    report["das_omp_rows"] = len(controlled)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true",
                        help="Maintainers only: record the current package files.")
    parser.add_argument("--shipwreck", action="store_true",
                        help="Also verify the extracted companion shipwreck asset.")
    parser.add_argument("--desktop-timing", action="store_true",
                        help="Also verify the extracted desktop MATLAB timing inputs.")
    args = parser.parse_args()
    report = science()
    computed = {p.relative_to(ROOT).as_posix(): digest(p) for p in files()}
    if args.write_manifest:
        MANIFEST.write_text("".join(f"{hash_}  {name}\n" for name, hash_ in computed.items()), encoding="utf-8")
    recorded = dict((name, hash_) for hash_, name in (
        line.split("  ", 1) for line in MANIFEST.read_text(encoding="utf-8").splitlines()))
    require(computed == recorded, "Manifest mismatch: " + ", ".join(sorted(
        k for k in computed.keys() | recorded.keys() if computed.get(k) != recorded.get(k))))
    if args.shipwreck:
        report["shipwreck_asset_files_verified"] = verify_asset(
            GEOMETRY, GEOMETRY_MANIFEST, "shipwreck-data_v10.zip")
    if args.desktop_timing:
        report["desktop_timing_asset_files_verified"] = verify_asset(
            DESKTOP, DESKTOP_MANIFEST, "desktop-timing-data_v10.zip")
    print(json.dumps({"status": "pass", "files_verified": len(computed), **report}, indent=2))


if __name__ == "__main__":
    main()
