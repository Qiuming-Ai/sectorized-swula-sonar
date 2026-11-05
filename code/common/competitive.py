"""Controlled physical wideband DAS/OMP comparison; isolated outputs only.

The principal algorithm comparison is paired: one broad LFM transmitter,
the same 32 physical random receivers, common Tx/Rx FE, and exactly the same
noisy measurements. Native CDM/FDM/Dense systems are a separate reference,
not a one-factor architecture or algorithm superiority experiment.
"""
from __future__ import annotations

import os
for _thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import argparse
from contextlib import nullcontext
import ctypes
import gc
import hashlib
import json
import platform
import shutil
import sys
import time

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
try:
    from threadpoolctl import threadpool_info, threadpool_limits
except ImportError:
    def threadpool_info():
        return [{"runtime_inspection": "threadpoolctl unavailable; fresh process sets BLAS/OpenMP environment to 1 before NumPy import"}]

    def threadpool_limits(limits=1):
        return nullcontext()

import numerical_model as model

ROOT = Path(__file__).resolve().parents[1]
NF = 8192
N_FREQUENCIES = 41
GRIDS = ((0.01, 0.5), (0.005, 0.25))
DISCREPANCY_FACTOR = 1.3
ITERATION_CAP = 32
NOISE_CALIBRATIONS = 30
SEEDS = (2026106300, 2026106301, 2026106302)


def build_receive_layout(cfg, min_spacing_width=False):
    """Same seed, first admissible candidates; never select by performance.

    Width-constrained mode repairs the formerly ideal, overlapping geometry.
    Both aperture endpoints are present before candidate rejection starts.
    Every candidate is a scalar draw from the original uniform distribution.
    """
    seed = 2026106001
    span = float(np.ptp(model.positions(cfg, "fdm")))
    width = float(cfg.element_width_lambda * cfg.c / cfg.fc)
    endpoints = [-span / 2, span / 2]
    rng = np.random.default_rng(seed)
    candidates = []
    if min_spacing_width:
        if span < (cfg.n_sparse - 1) * width:
            raise ValueError("Declared aperture cannot fit the fixed element-width constraint")
        accepted = list(endpoints)
        for attempt in range(1, 1000001):
            candidate = float(rng.uniform(-span / 2, span / 2))
            minimum = float(np.min(np.abs(np.asarray(accepted) - candidate)))
            admissible = minimum >= width
            candidates.append({"uniform_candidate_number": attempt, "position_m": candidate,
                               "minimum_distance_to_previous_accepted_m": minimum,
                               "accepted": bool(admissible)})
            if admissible:
                accepted.append(candidate)
            if len(accepted) == cfg.n_sparse:
                break
        else:
            raise RuntimeError("Fixed deterministic rejection sampling exhausted its safety cap")
        positions = np.sort(np.asarray(accepted))
    else:
        draw = rng.uniform(-span / 2, span / 2, cfg.n_sparse - 2)
        positions = np.sort(np.r_[endpoints, draw])
        candidates = [{"uniform_candidate_number": i + 1, "position_m": float(value),
                       "minimum_distance_to_previous_accepted_m": None, "accepted": True}
                      for i, value in enumerate(draw)]
    spacing = np.diff(positions)
    layout = {"layout_mode": "fixed_seed_element_width_rejection" if min_spacing_width else "original_ideal_unrestricted_sampling",
              "collision_correction_not_performance_tuning": bool(min_spacing_width),
              "fixed_seed": seed, "receive_channels": len(positions), "aperture_span_m": span,
              "endpoint_positions_m": endpoints, "element_width_lambda": cfg.element_width_lambda,
              "element_width_m": width, "minimum_required_spacing_m": width if min_spacing_width else None,
              "minimum_actual_center_spacing_m": float(spacing.min()),
              "adjacent_gaps_less_than_element_width": int(np.count_nonzero(spacing < width)),
              "attempted_uniform_candidates": len(candidates),
              "accepted_uniform_candidates": cfg.n_sparse - 2,
              "rejected_uniform_candidates": len(candidates) - (cfg.n_sparse - 2),
              "rule": "Endpoints fixed first; sequential uniform candidates rejected iff distance to ANY accepted center is <0.45lambda; accept first 30 admissible candidates; sort only after selection" if min_spacing_width else "original first 30 unrestricted uniform candidates plus two endpoints",
              "positions_m": positions.tolist(), "positions_sha256": digest_array(positions),
              "candidate_sequence": candidates,
              "manufacturing_scope": "minimum-width center packing only; no mutual-coupling, mounting, fabrication, or hardware validation"}
    if min_spacing_width and layout["adjacent_gaps_less_than_element_width"]:
        raise AssertionError("Width-constrained layout still has overlap")
    return positions, layout


def memory_state():
    """Actual process memory, not estimated dictionary bytes or GPU memory."""
    result = {"process_working_set_bytes": None,
              "process_peak_working_set_bytes": None,
              "memory_measurement_scope": "process lifetime; includes imports and earlier tasks"}
    if sys.platform == "win32":
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                (name, ctypes.c_size_t) for name in (
                    "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                    "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                    "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]

        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        if psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            result["process_working_set_bytes"] = int(counters.WorkingSetSize)
            result["process_peak_working_set_bytes"] = int(counters.PeakWorkingSetSize)
    else:
        try:
            import psutil
            result["process_working_set_bytes"] = int(psutil.Process().memory_info().rss)
            import resource
            scale = 1 if sys.platform == "darwin" else 1024
            result["process_peak_working_set_bytes"] = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale)
        except ImportError:
            pass
    return result


def digest_array(array):
    return hashlib.sha256(np.ascontiguousarray(array).view(np.uint8)).hexdigest()


def grid(dr, da):
    rr = np.arange(30.3, 30.700001, dr)
    aa = np.arange(-50.0, 50.000001, da)
    rg, ag = np.meshgrid(rr, aa, indexing="ij")
    return rr, aa, rg.ravel(), ag.ravel()


def physical_atoms(cfg, xpos, frequencies, effective_lfm_spectrum, ranges, angles):
    """Exact spherical paths, full RF f=fc+nu, actual Rx element pattern.

    The broad ideal point Tx has unit directivity; effective_lfm_spectrum
    contains the actual, common-FE-shaped Tx and common Rx FE transfer.
    """
    r = np.atleast_1d(ranges)
    a = np.deg2rad(np.atleast_1d(angles))
    x = r * np.sin(a)
    z = r * np.cos(a)
    distance = np.hypot(x[None, :] - xpos[:, None], z[None, :])
    arrival_u = (x[None, :] - xpos[:, None]) / distance
    delay = (r[None, :] + distance - 2 * cfg.gate_range) / cfg.c
    directivity = np.sinc(cfg.element_width_lambda / cfg.fc * frequencies[:, None, None] * arrival_u[None, :, :])
    atoms = (effective_lfm_spectrum[:, None, None] * directivity
             * np.exp(-2j * np.pi * frequencies[:, None, None] * delay[None, :, :]))
    return atoms.reshape(len(frequencies) * len(xpos), len(r))


def build_dictionary(cfg, xpos, frequencies, effective_spectrum, ranges, angles):
    start = time.perf_counter()
    cpu_start = time.process_time()
    matrix = np.empty((len(frequencies) * len(xpos), len(ranges)), dtype=np.complex128)
    for left in range(0, len(ranges), 256):
        right = min(left + 256, len(ranges))
        matrix[:, left:right] = physical_atoms(cfg, xpos, frequencies, effective_spectrum,
                                              ranges[left:right], angles[left:right])
    norms = np.linalg.norm(matrix, axis=0)
    if np.any(norms <= 0):
        raise ValueError("Zero-energy physical atom")
    matrix /= norms[None, :]
    adjoint = np.ascontiguousarray(matrix.conj().T)
    report = {"dictionary_build_seconds_local_CPU_wall": time.perf_counter() - start,
              "dictionary_build_seconds_process_CPU": time.process_time() - cpu_start,
              "dictionary_bytes": matrix.nbytes,
              "conjugate_dictionary_bytes": adjoint.nbytes,
              "dictionary_plus_conjugate_bytes": matrix.nbytes + adjoint.nbytes,
              "dictionary_atoms": len(ranges), "dictionary_dtype": str(matrix.dtype),
              **memory_state()}
    return matrix, adjoint, norms, report


def omp(matrix, adjoint, norms, measurement, sigma):
    """K is never an input. Initial discrepancy check permits empty scenes."""
    started = time.perf_counter()
    cpu_started = time.process_time()
    residual = measurement.copy()
    selected = []
    coefficients = np.empty(0, dtype=complex)
    tolerance = DISCREPANCY_FACTOR * sigma * np.sqrt(len(measurement))
    trace = [{"iteration": 0, "selected_dictionary_index": None,
              "residual_norm": float(np.linalg.norm(residual)),
              "stop_tolerance": float(tolerance),
              "residual_over_discrepancy_threshold": float(np.linalg.norm(residual) / tolerance)}]
    while np.linalg.norm(residual) > tolerance and len(selected) < ITERATION_CAP:
        score = np.abs(adjoint @ residual)
        if selected:
            score[selected] = -np.inf
        selected.append(int(np.argmax(score)))
        coefficients = np.linalg.lstsq(matrix[:, selected], measurement, rcond=None)[0]
        residual = measurement - matrix[:, selected] @ coefficients
        trace.append({"iteration": len(selected), "selected_dictionary_index": selected[-1],
                      "residual_norm": float(np.linalg.norm(residual)),
                      "stop_tolerance": float(tolerance),
                      "residual_over_discrepancy_threshold": float(np.linalg.norm(residual) / tolerance)})
    image = np.zeros(len(norms), dtype=float)
    if selected:
        image[selected] = np.abs(coefficients / norms[selected])
    converged = bool(np.linalg.norm(residual) <= tolerance)
    report = {"solver_seconds_local_CPU_wall": time.perf_counter() - started,
              "solver_seconds_process_CPU": time.process_time() - cpu_started,
              "OMP_iterations": len(selected), "discrepancy_stop_met": converged,
              "stop_reason": "discrepancy" if converged else "maxiter",
              "cap_reached_without_discrepancy_convergence": bool(not converged and len(selected) == ITERATION_CAP),
              "residual_over_discrepancy_threshold": float(np.linalg.norm(residual) / tolerance),
              **memory_state()}
    return image, report, trace


def das(adjoint, norms, measurement):
    """Normalized physical-model matched DAS; coherent sum over all samples."""
    start = time.perf_counter()
    cpu_start = time.process_time()
    image = np.abs(adjoint @ measurement) / norms
    return image, {"solver_seconds_local_CPU_wall": time.perf_counter() - start,
                   "solver_seconds_process_CPU": time.process_time() - cpu_start,
                   "OMP_iterations": None, "discrepancy_stop_met": None,
                   "stop_reason": "not_iterative", "cap_reached_without_discrepancy_convergence": None,
                   "residual_over_discrepancy_threshold": None, **memory_state()}


def blind_detection(image, rr, aa, threshold, truth_r, truth_a):
    """Global image candidates first; truth only enters final assignment."""
    amplitude = np.abs(image)
    nr = max(1, int(np.ceil(0.03 / (rr[1] - rr[0]))))
    na = max(1, int(np.ceil(2.0 / (aa[1] - aa[0]))))
    local = ndimage.maximum_filter(amplitude, size=(2 * nr + 1, 2 * na + 1),
                                  mode="constant", cval=-np.inf)
    indices = np.argwhere((amplitude == local) & (amplitude >= threshold) & (amplitude > 0))
    order = np.argsort(amplitude[indices[:, 0], indices[:, 1]])[::-1] if len(indices) else []
    candidates = []
    for i, j in indices[order] if len(indices) else []:
        if any(abs(rr[i] - r) <= 0.03 and abs(aa[j] - a) <= 2.0 for r, a, _ in candidates):
            continue
        candidates.append((float(rr[i]), float(aa[j]), float(amplitude[i, j])))
    matches = []
    if candidates and len(truth_r):
        cr = np.array([row[0] for row in candidates])
        ca = np.array([row[1] for row in candidates])
        rd = np.abs(cr[:, None] - truth_r[None, :])
        ad = np.abs(ca[:, None] - truth_a[None, :])
        valid = (rd <= 0.02) & (ad <= 1.0)
        cost = np.where(valid, (rd / 0.02) ** 2 + ad ** 2, 1e6)
        pi, ti = linear_sum_assignment(cost)
        matches = [(int(p), int(t)) for p, t in zip(pi, ti) if valid[p, t]]
    return {"detection_probability": float(len(matches) / len(truth_r)) if len(truth_r) else None,
            "true_detections": len(matches), "false_peaks": len(candidates) - len(matches),
            "blind_peak_count": len(candidates),
            "matched_range_error_rmse_m": float(np.sqrt(np.mean([(candidates[p][0] - truth_r[t]) ** 2 for p, t in matches]))) if matches else None,
            "matched_angle_error_rmse_deg": float(np.sqrt(np.mean([(candidates[p][1] - truth_a[t]) ** 2 for p, t in matches]))) if matches else None}, candidates


def calibrated_threshold(adjoint, norms, sigma, seed):
    """Frozen common threshold from independent noise-only DAS maxima.

    Shared by paired DAS/OMP; not a claim of equal calibrated OMP false-alarm
    probability. OMP discrepancy stopping rejects ordinary noise-only inputs.
    """
    maxima = []
    for repetition in range(NOISE_CALIBRATIONS):
        rng = np.random.default_rng(seed + repetition)
        noise = sigma * (rng.standard_normal(adjoint.shape[1]) + 1j * rng.standard_normal(adjoint.shape[1])) / np.sqrt(2)
        maxima.append(float(np.max(np.abs(adjoint @ noise) / norms)))
    return float(np.quantile(maxima, 0.95)), maxima


def streamed_das(cfg, xpos, frequencies, spectrum, ranges, angles, measurement):
    """Geometry-only control without retaining a second enormous dictionary."""
    start = time.perf_counter()
    result = np.empty(len(ranges), dtype=float)
    for left in range(0, len(ranges), 256):
        right = min(left + 256, len(ranges))
        atoms = physical_atoms(cfg, xpos, frequencies, spectrum, ranges[left:right], angles[left:right])
        result[left:right] = np.abs(atoms.conj().T @ measurement) / np.sum(np.abs(atoms) ** 2, axis=0)
    return result, time.perf_counter() - start


def native_system_reference(cfg, run, truth_r, truth_a, coefficients, sigma, bins, xpos):
    """Noisy native references use pre-FE frequency noise; separately scoped."""
    rr, aa, _, _ = grid(*GRIDS[0])
    rows, arrays = [], {}
    send, matched, nu, _ = model.spectra(cfg, "dense", NF)
    f = cfg.fc + nu
    random_raw = physical_atoms(cfg, xpos, f, send[:, 0], truth_r, truth_a) @ coefficients
    random_raw = random_raw.reshape(NF, len(xpos))
    # Include a random-Rx LFM reference through identical continuous Rx FE.
    # This general custom-position implementation copies the mathematical
    # core DAS, not its native geometry selector; no old code is modified.
    def custom_reconstruction(raw, positions):
        u = np.sin(np.deg2rad(aa))
        mf = np.fft.ifft(raw * model.frontend_response(cfg, NF)[:, None] * np.conj(matched[:, 0, None]), axis=0)
        fine = model.signal.resample(mf, NF * cfg.interpolation_oversampling, axis=0)
        image = np.zeros((len(rr), len(aa)), dtype=complex)
        weights = model.weights(len(positions))
        xx = rr[:, None] * u[None, :]
        zz = rr[:, None] * np.cos(np.deg2rad(aa))[None, :]
        for j, position in enumerate(positions):
            tau = (rr[:, None] + np.hypot(xx - position, zz)) / cfg.c - 2 * cfg.gate_range / cfg.c
            index = tau * cfg.fs * cfg.interpolation_oversampling
            low = np.floor(index).astype(int)
            fraction = index - low
            values = fine[low, j] * (1 - fraction) + fine[low + 1, j] * fraction
            image += weights[j] * values * np.exp(2j * np.pi * cfg.fc * tau)
        return image

    for method in ("random_lfm", "dense", "cdm", "fdm"):
        started = time.perf_counter()
        if method == "random_lfm":
            raw, mf, metadata = random_raw, matched, {"method": "broad LFM / 32 random Rx"}
        else:
            raw, mf, metadata = model.point_raw(cfg, method, truth_r, truth_a, coefficients,
                                              geometry="spherical", nfft=NF)
        synthesis = time.perf_counter() - started
        for repetition, seed in enumerate(SEEDS):
            rng = np.random.default_rng(seed)
            standardized = (rng.standard_normal((NF, 128)) + 1j * rng.standard_normal((NF, 128))) / np.sqrt(2)
            # Independent frequency-domain CN noise is exactly equivalent to
            # time-domain white noise with std sigma/sqrt(NF), before Rx FE.
            noisy = raw + sigma * standardized[:, :raw.shape[1]]
            started = time.perf_counter()
            image = (custom_reconstruction(noisy, xpos) if method == "random_lfm"
                     else model.reconstruct_points(cfg, method, noisy, mf, rr, aa,
                                                   focus_model="spherical"))
            solve = time.perf_counter() - started
            metrics, _ = blind_detection(image, rr, aa, 0.05 * np.max(np.abs(image)), truth_r, truth_a)
            rows.append({"comparison_family": "noisy_native_system_reference_not_single_factor",
                         "method": method, "noise_replication": repetition, "noise_seed": seed,
                         "receive_channels": raw.shape[1], "preFE_frequency_noise_sigma": sigma,
                         "equivalent_preFE_time_noise_sigma": sigma / np.sqrt(NF),
                         "selected_frequency_input_snr_db": float(10 * np.log10(np.mean(np.abs(raw[bins]) ** 2) / sigma ** 2)),
                         "selected_frequency_postFE_snr_db": float(10 * np.log10(np.mean(np.abs(raw[bins] * model.frontend_response(cfg, NF)[bins, None]) ** 2)
                                                                                / np.mean(np.abs(model.frontend_response(cfg, NF)[bins]) ** 2) / sigma ** 2)),
                         "synthesis_seconds_local_CPU_wall": synthesis,
                         "solver_seconds_local_CPU_wall": solve,
                         "total_synthesis_plus_solver_seconds_local_CPU_wall": synthesis + solve,
                         "threshold_kind": "5pct_image_peak_diagnostic_not_equal_false_alarm",
                         **metrics, **memory_state()})
            arrays[f"{method}_rep{repetition}"] = np.abs(image).astype(np.float32)
            print("native noisy", method, repetition, metrics, flush=True)
        model.dump(run / f"native_{method}_metadata.json", metadata)
    model.write_csv(run / "native_system_metrics.csv", rows)
    np.savez_compressed(run / "native_system_images.npz", ranges_m=rr, angles_deg=aa, **arrays)
    return rows


def run_experiment(min_spacing_width=False):
    if not (ROOT / "core_ready.json").exists():
        raise RuntimeError("Core readiness marker missing; run the fixed model gates first")
    readiness = json.loads((ROOT / "core_ready.json").read_text(encoding="utf-8"))
    if readiness.get("gate") != "passed" or readiness.get("model_sha256") != model.sha(ROOT / "common" / "numerical_model.py"):
        raise RuntimeError("Core readiness hash does not match current numerical source")
    parent = ROOT.parent / "outputs" / "competitive"
    run = parent / ("run_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
    run.mkdir(parents=True, exist_ok=False)
    cfg = model.Config()
    started = time.perf_counter()
    random_positions, layout = build_receive_layout(cfg, min_spacing_width=min_spacing_width)
    model.dump(run / "receive_layout_frozen_before_run.json", layout)
    send, _, nu, tx_metadata = model.spectra(cfg, "dense", NF)
    bins = np.sort(np.random.default_rng(2026106002).choice(np.where(np.abs(nu) <= 45000)[0], N_FREQUENCIES, replace=False))
    frequencies = cfg.fc + nu[bins]
    spectrum = send[bins, 0] * model.frontend_response(cfg, NF)[bins]
    protocol = {"status": "frozen_before_training_and_test",
                "configuration": asdict(cfg), "grids": [{"range_step_m": dr, "angle_step_deg": da} for dr, da in GRIDS],
                "dictionary_domain": {"range_m": [30.3, 30.7], "angle_deg": [-50, 50]},
                "OMP_stop_rule": {"initial_discrepancy_check": True, "factor": DISCREPANCY_FACTOR, "maxiter": ITERATION_CAP},
                "source_count_given_to_solver": False,
                "training": "separate two-target scene; one sigma from its clean random-Rx measurement defines 20 dB; frozen before all test data",
                "test": "five continuous off-grid targets; three paired independent noise seeds; no oracle neighborhood search",
                "noise_seed_list": list(SEEDS), "array_seed": 2026106001, "frequency_seed": 2026106002,
                "selected_frequency_bins": bins.tolist(), "selected_baseband_hz": nu[bins].tolist(),
                "random_positions_m": random_positions.tolist(),
                "receive_layout": layout,
                "layout_correction_scope": "repair physical center overlap with fixed-seed deterministic rejection; no waveform, dictionary grid, noise seed, stopping, or performance-driven geometry tuning",
                "core_ready_marker": readiness,
                "controlled_noise": "additive white CN(0,sigma^2) on selected POST-FE frequency-channel measurements; not a preFE full-record noise experiment",
                "native_reference_noise": "same numerical sigma in all preFE FFT bins; equivalent white full-record sample sigma/sqrt(8192); Rx FE then filters noise",
                "native_reference_not_controlled_ranking": True,
                "controlled_DAS": "normalized physical-model adjoint coherent DAS over identical 41x32 measurements",
                "dictionary": "exact spherical Tx point + physical Rx path; full RF frequencies; Rx element sinc pattern; actual Tx FE and actual Rx FE",
                "detection": {"candidates": "global 2D NMS, no truth access", "NMS_range_m": 0.03, "NMS_angle_deg": 2,
                              "matching": "one-to-one minimum-cost assignment after candidates", "range_tolerance_m": 0.02, "angle_tolerance_deg": 1,
                              "threshold": "95th percentile of 30 independent noise-only global DAS amplitude maxima; same absolute threshold for paired DAS and OMP; NOT equal OMP false-alarm calibration",
                              "secondary": "5pct-of-own-peak diagnostic retained separately"},
                "timing": "measured local CPU wall and process time; no Jetson, GPU, or end-to-end hardware throughput claim",
                "memory": "dictionary and conjugate exact bytes AND measured lifetime peak process working set; neither is GPU VRAM",
                "software": {"python": sys.version, "platform": platform.platform(), "threadpools": threadpool_info()},
                "source_sha256_before": {"competitive.py": model.sha(Path(__file__)), "numerical_model.py": model.sha(ROOT / "common" / "numerical_model.py")}}
    model.dump(run / "protocol_before_run.json", protocol)
    # Independent training must precede test construction and evaluation.
    train_r, train_a = np.array([30.36, 30.64]), np.array([-18.0, 27.0])
    train_c = np.exp(1j * np.random.default_rng(2026106100).uniform(-np.pi, np.pi, 2)) / np.sqrt(2)
    training_clean = physical_atoms(cfg, random_positions, frequencies, spectrum, train_r, train_a) @ train_c
    sigma = float(np.linalg.norm(training_clean) / np.sqrt(training_clean.size) / 10)
    model.dump(run / "training_sigma_frozen_before_test.json", {"sigma": sigma, "training_snr_db": 20.0,
               "training_truth_ranges_m": train_r.tolist(), "training_truth_angles_deg": train_a.tolist(),
               "training_clean_sha256": digest_array(training_clean), "sigma_reused_unchanged_for_test_and_all_geometries": True})
    test_r = np.array([30.3537, 30.4213, 30.4967, 30.5731, 30.6459])
    test_a = np.array([-39.37, -17.23, 0.37, 24.13, 46.41])
    test_c = np.exp(1j * np.random.default_rng(2026106101).uniform(-np.pi, np.pi, 5)) / np.sqrt(5)
    synthesis_start = time.perf_counter()
    test_clean = physical_atoms(cfg, random_positions, frequencies, spectrum, test_r, test_a) @ test_c
    synthesis_time = time.perf_counter() - synthesis_start
    measurements = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        noise = (rng.standard_normal((N_FREQUENCIES, 128)) + 1j * rng.standard_normal((N_FREQUENCIES, 128))) / np.sqrt(2)
        measurements.append((test_clean.reshape(N_FREQUENCIES, 32) + sigma * noise[:, :32]).ravel())
    # Validate the custom physical generator against the actual core native
    # Dense point-Rx model before interpreting the competitive reconstructions.
    core_raw, _, _ = model.point_raw(cfg, "dense", test_r, test_a, test_c, geometry="spherical", nfft=NF)
    dense_positions = model.positions(cfg, "dense")
    custom_dense = physical_atoms(cfg, dense_positions, frequencies, spectrum, test_r, test_a) @ test_c
    expected = (core_raw[bins] * model.frontend_response(cfg, NF)[bins, None]).ravel()
    relative_error = float(np.linalg.norm(custom_dense - expected) / np.linalg.norm(expected))
    model.dump(run / "custom_model_gate.json", {"relative_error": relative_error, "tolerance": 1e-10,
               "passed": relative_error <= 1e-10, "comparison": "custom spherical physical selected-FE data vs core Dense raw FFT through common Rx FE"})
    if relative_error > 1e-10:
        raise RuntimeError("Custom physical dictionary/data model disagrees with core")
    del core_raw, custom_dense, expected
    rows, traces, peaks, calibration_rows = [], [], [], []
    saved = {"test_ranges_m": test_r, "test_angles_deg": test_a,
             "random_positions_m": random_positions, "selected_baseband_hz": nu[bins],
             "controlled_clean_measurements": test_clean, "controlled_noisy_measurements": np.stack(measurements)}
    test_snr = float(10 * np.log10(np.mean(np.abs(test_clean) ** 2) / sigma ** 2))
    for grid_index, (dr, da) in enumerate(GRIDS):
        rr, aa, gr, ga = grid(dr, da)
        matrix, adjoint, norms, build_report = build_dictionary(cfg, random_positions, frequencies, spectrum, gr, ga)
        threshold, maxima = calibrated_threshold(adjoint, norms, sigma, 2026106400 + grid_index * 100)
        calibration_rows.extend({"range_step_m": dr, "angle_step_deg": da, "repetition": i,
                                 "noise_only_DAS_maximum_amplitude": value} for i, value in enumerate(maxima))
        model.dump(run / f"calibration_frozen_{grid_index}.json", {"threshold": threshold,
                   "sigma": sigma, "noise_trials": NOISE_CALIBRATIONS, "no_test_data_used_for_threshold": True,
                   "same_threshold_applied_to_DAS_and_OMP": True, "equal_false_alarm_guarantee": False})
        print("dictionary ready", dr, da, build_report, "threshold", threshold, flush=True)
        # Training outcome is retained but never used to adjust stopping rules.
        rng = np.random.default_rng(2026106200)
        train_y = training_clean + sigma * (rng.standard_normal(training_clean.size) + 1j * rng.standard_normal(training_clean.size)) / np.sqrt(2)
        train_im, train_report, train_trace = omp(matrix, adjoint, norms, train_y, sigma)
        model.dump(run / f"training_OMP_{grid_index}.json", train_report)
        traces.extend({"scene": "independent_training_2", "range_step_m": dr, "angle_step_deg": da,
                       "noise_replication": 0, **entry} for entry in train_trace)
        for repetition, measurement in enumerate(measurements):
            for algorithm in ("DAS", "OMP"):
                if algorithm == "DAS":
                    image, solve_report = das(adjoint, norms, measurement)
                    trace = []
                else:
                    image, solve_report, trace = omp(matrix, adjoint, norms, measurement, sigma)
                image = image.reshape(len(rr), len(aa))
                metrics, candidates = blind_detection(image, rr, aa, threshold, test_r, test_a)
                diagnostic, _ = blind_detection(image, rr, aa, 0.05 * np.max(image), test_r, test_a)
                row = {"comparison_family": "paired_algorithm_same_Tx_Rx_FE_noisy_data",
                       "scene": "offgrid_test_5", "method": algorithm, "range_step_m": dr, "angle_step_deg": da,
                       "noise_replication": repetition, "noise_seed": SEEDS[repetition], "sigma_frequency_samples_postFE": sigma,
                       "actual_test_selected_frequency_snr_db": test_snr,
                       "measurement_sha256": digest_array(measurement), "absolute_detection_threshold": threshold,
                       "threshold_scope": "common noise-DAS-global-max calibration; not equal OMP Pfa",
                       "synthesis_seconds_local_CPU_wall": synthesis_time, **build_report, **solve_report,
                       "total_unamortized_synthesis_build_solve_seconds_local_CPU_wall": synthesis_time + build_report["dictionary_build_seconds_local_CPU_wall"] + solve_report["solver_seconds_local_CPU_wall"],
                       "total_scope": "numerical synthesis + one dictionary build + one solve; excludes noise calibration, metrics, serialization",
                       **metrics, **{"relative5pct_diagnostic_" + key: value for key, value in diagnostic.items()}}
                rows.append(row)
                traces.extend({"scene": "offgrid_test_5", "range_step_m": dr, "angle_step_deg": da,
                               "noise_replication": repetition, **entry} for entry in trace)
                peaks.extend({"method": algorithm, "range_step_m": dr, "angle_step_deg": da,
                              "noise_replication": repetition, "candidate_range_m": r,
                              "candidate_angle_deg": a, "candidate_amplitude": value} for r, a, value in candidates)
                saved[f"{algorithm}_grid{grid_index}_rep{repetition}"] = image.astype(np.float32)
                saved[f"ranges_grid{grid_index}"] = rr
                saved[f"angles_grid{grid_index}"] = aa
                print("controlled", dr, da, repetition, algorithm, solve_report["stop_reason"], metrics, flush=True)
        del matrix, adjoint, norms, train_im
        gc.collect()
    model.write_csv(run / "controlled_metrics.csv", rows)
    model.write_csv(run / "iteration_traces.csv", traces)
    model.write_csv(run / "blind_candidates.csv", peaks)
    model.write_csv(run / "noise_threshold_calibration.csv", calibration_rows)
    # Broad-LFM receive-geometry control: all three arrays share transmitter,
    # common FE and the same sigma. Only Rx geometry/count varies.
    geometry_rows = []
    rr, aa, gr, ga = grid(*GRIDS[0])
    for name, xpos in (("32random", random_positions), ("32SWULA", model.positions(cfg, "fdm")),
                       ("128Dense", model.positions(cfg, "dense"))):
        clean = physical_atoms(cfg, xpos, frequencies, spectrum, test_r, test_a) @ test_c
        for repetition, seed in enumerate(SEEDS):
            rng = np.random.default_rng(seed)
            noise = (rng.standard_normal((N_FREQUENCIES, 128)) + 1j * rng.standard_normal((N_FREQUENCIES, 128))) / np.sqrt(2)
            measurement = (clean.reshape(N_FREQUENCIES, len(xpos)) + sigma * noise[:, :len(xpos)]).ravel()
            image, elapsed = streamed_das(cfg, xpos, frequencies, spectrum, gr, ga, measurement)
            image = image.reshape(len(rr), len(aa))
            metrics, _ = blind_detection(image, rr, aa, 0.05 * image.max(), test_r, test_a)
            geometry_rows.append({"comparison_family": "geometry_control_same_broadLFM_FE_sigma_DAS",
                                  "receive_geometry": name, "receive_channels": len(xpos), "noise_replication": repetition,
                                  "sigma_frequency_samples_postFE": sigma,
                                  "selected_frequency_snr_db": float(10 * np.log10(np.mean(np.abs(clean) ** 2) / sigma ** 2)),
                                  "dictionary_build_plus_solver_seconds_local_CPU_wall": elapsed,
                                  "threshold_kind": "5pct_image_peak_diagnostic_not_equal_false_alarm", **metrics})
            saved[f"geometry_{name}_rep{repetition}"] = image.astype(np.float32)
        print("geometry controlled", name, flush=True)
    model.write_csv(run / "geometry_control_metrics.csv", geometry_rows)
    native_rows = native_system_reference(cfg, run, test_r, test_a, test_c, sigma, bins, random_positions)
    np.savez_compressed(run / "controlled_raw_images_and_measurements.npz", **saved)
    model.write_csv(run / "truth_catalog.csv", [{"range_m": float(r), "angle_deg": float(a),
                    "amplitude": float(abs(c)), "phase_rad": float(np.angle(c))}
                    for r, a, c in zip(test_r, test_a, test_c)])
    model.dump(run / "scope.json", {"primary_inference": "paired DAS versus OMP on identical broadLFM/32random/commonFE/postFE-noisy measurements",
               "receive_layout": {key: value for key, value in layout.items() if key != "candidate_sequence"},
               "secondary_geometry_control": "same broad point LFM and FE; 32random vs 32SWULA vs 128Dense; channel count also changes for Dense",
               "native_system_reference": "32random-LFM vs 128Dense-LFM vs 32SWULA-CDM/FDM; same preFE noise sigma; different Tx and geometry; not a single-factor or algorithm superiority ranking",
               "source_count_unknown": True, "training_count_2_test_count_5": True,
               "sigma_training_locked": sigma, "three_noise_repetitions_one_fixed_test_scene": True,
               "no_truth_centered_candidate_search": True, "coarse_and_fine_results_all_retained": True,
               "native_and_geometry_relative_peak_metrics": "diagnostic only; not calibrated sensitivity at equal false alarms",
               "OMP_threshold": "shared DAS-calibrated absolute threshold, with false peaks reported; no claim of equal Pfa",
               "memory_scope": "real process working-set lifetime peak includes imports/previous phases; dictionary bytes alone underestimate working memory",
               "CPU_scope": "all measured local CPU; numerical stages only, not Jetson/GPU or acquired-frame end-to-end performance",
               "tx_metadata": tx_metadata})
    snapshot = run / "code_snapshot"
    snapshot.mkdir()
    for path in (Path(__file__), ROOT / "common" / "numerical_model.py"):
        shutil.copy2(path, snapshot / path.name)
    source_after = {"competitive.py": model.sha(Path(__file__)), "numerical_model.py": model.sha(ROOT / "common" / "numerical_model.py")}
    if source_after != protocol["source_sha256_before"]:
        raise RuntimeError("Source changed during execution; results not accepted")
    fine_omp = [row for row in rows if row["method"] == "OMP" and row["range_step_m"] == GRIDS[1][0]]
    coarse_omp = [row for row in rows if row["method"] == "OMP" and row["range_step_m"] == GRIDS[0][0]]
    receipt = {"status": "computed_pending_author_validation", "source_sha256": source_after,
               "layout_mode": layout["layout_mode"],
               "collision_correction_not_performance_tuning": bool(min_spacing_width),
               "minimum_actual_center_spacing_m": layout["minimum_actual_center_spacing_m"],
               "scope": "simulation only; no superiority claim inferred automatically",
               "fine_grid_all_discrepancy_stops_met": all(row["discrepancy_stop_met"] for row in fine_omp),
               "coarse_grid_all_discrepancy_stops_met": all(row["discrepancy_stop_met"] for row in coarse_omp),
               "controlled_metric_rows": len(rows), "geometry_metric_rows": len(geometry_rows), "native_metric_rows": len(native_rows),
               "total_experiment_seconds_local_CPU_wall": time.perf_counter() - started, **memory_state(),
               "artifacts": [{"path": str(path.relative_to(run)), "sha256": model.sha(path)}
                             for path in sorted(run.rglob("*")) if path.is_file()]}
    model.dump(run / "receipt.json", receipt)
    model.dump(parent / "latest.json", {"path": str(run.relative_to(ROOT.parent)), "status": receipt["status"]})
    if min_spacing_width:
        model.dump(parent / "primary_minspacing_latest.json", {"path": str(run.relative_to(ROOT.parent)),
                   "status": receipt["status"], "layout_mode": layout["layout_mode"],
                   "reason": "physical center-overlap correction, not selection for favorable performance"})
    print(json.dumps({"run": str(run), **{key: value for key, value in receipt.items() if key != "artifacts"}}, ensure_ascii=False), flush=True)
    return run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-spacing-width", action="store_true",
                        help="Fixed-seed rejection layout with minimum spacing equal to original 0.45lambda element width")
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        run_experiment(min_spacing_width=args.min_spacing_width)


if __name__ == "__main__":
    main()
