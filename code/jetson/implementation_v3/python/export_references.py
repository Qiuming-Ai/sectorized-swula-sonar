"""Export frozen candidate banks and FP64 CPU reference images.

Use --fixtures to produce a development subset, then --resume for the unchanged
remaining cases. Existing numerical outputs are not overwritten. Every image is
also compared with the frozen v2 full-rate FP64 reference, without a gain fit.
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import numpy as np
from fdm_reference import (HERE, V2_PYTHON, CANDIDATES, FIXTURE_NAMES, LUT_PHASES,
    KAISER_BETA, candidate_bank, ensure_applicable, image_from_mf, read_method,
    read_raw, receive_fft, matched_time, comparison, save_binary, sha, dump)


def export_static(root, bank, name):
    candidate = candidate_bank(bank, name)
    dest = root / name
    dest.mkdir(parents=True, exist_ok=True)
    meta = dict(schema_version=3, **CANDIDATES[name],
                lut_phases=LUT_PHASES, kaiser_beta=KAISER_BETA,
                low_samples=candidate["low_samples"], fs_low=candidate["fs_low"],
                fir_delay_compensated_samples=candidate["fir_delay_compensated_samples"],
                native_fft_samples=int(bank["meta"]["fft_samples"]),
                shared_bank=str(V2_PYTHON / "generated" / "fdm" / "static"))
    content = "".join(f"{key}={value}\n" for key, value in meta.items())
    meta_path = dest / "meta.txt"
    if meta_path.exists() and meta_path.read_text(encoding="utf-8") != content:
        raise RuntimeError(f"Existing static candidate metadata differs: {name}")
    if not meta_path.exists():
        meta_path.write_text(content, encoding="utf-8")
    files = []
    for filename, arr, dtype in (("H.bin", candidate["h"], "<c8"),
            ("H64.bin", candidate["h"], "<c16"), ("mix.bin", candidate["mix"], "<c8"),
            ("mix64.bin", candidate["mix"], "<c16"), ("offsets.bin", candidate["offsets"], "<f8"),
            ("lut.bin", candidate["lut"], "<f4"), ("lut64.bin", candidate["lut"], "<f8")):
        path = dest / filename
        if path.exists():
            expected = np.ascontiguousarray(arr, dtype=dtype)
            actual = np.fromfile(path, dtype=dtype)
            if actual.shape != (expected.size,) or not np.array_equal(actual, expected.ravel()):
                raise RuntimeError(f"Refusing to overwrite changed static {path}")
            files.append(dict(path=str(path), dtype=np.dtype(dtype).str, shape=list(expected.shape), bytes=path.stat().st_size, sha256=sha(path)))
        else:
            files.append(save_binary(path, arr, dtype))
    return candidate, files


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "generated")
    parser.add_argument("--fixtures", nargs="+", choices=FIXTURE_NAMES, default=FIXTURE_NAMES)
    parser.add_argument("--candidates", nargs="+", choices=list(CANDIDATES), default=list(CANDIDATES))
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root = args.out.resolve()
    if root.exists() and any(root.iterdir()) and not args.resume:
        raise FileExistsError("Output already exists; --resume verifies unchanged static bank")
    root.mkdir(parents=True, exist_ok=True)
    shared = V2_PYTHON / "generated" / "fdm"
    bank = read_method(shared / "static")
    ensure_applicable(bank)
    candidates, static_files = {}, []
    for name in args.candidates:
        candidates[name], files = export_static(root, bank, name)
        static_files.extend(files)
    audit = []
    for fixture in args.fixtures:
        needed = [name for name in candidates if not (root / name / "cases" / fixture / "comparison_to_v2.json").exists()]
        if not needed:
            print(f"SKIP existing complete {fixture}", flush=True)
            continue
        raw_path = shared / fixture / "raw.bin"
        raw = read_raw(raw_path, bank["meta"])
        fullrate_path = shared / fixture / "reference.bin"
        fullrate = np.fromfile(fullrate_path, dtype="<c16").reshape(len(bank["ranges"]), len(bank["angles"]))
        fe = receive_fft(raw, bank)
        matched = [matched_time(fe, bank, k) for k in range(int(bank["meta"]["waveforms"]))]
        for name in needed:
            began = time.perf_counter()
            candidate = candidates[name]
            image = image_from_mf(matched, bank, candidate)
            dest = root / name / "cases" / fixture
            dest.mkdir(parents=True, exist_ok=True)
            artifact = save_binary(dest / "reference.bin", image, "<c16")
            stats = comparison(fullrate, image, bank)
            stats.update(candidate=name, fixture=fixture,
                         cpu_reference_seconds=time.perf_counter() - began,
                         reference_fullrate_path=str(fullrate_path), reference_fullrate_sha256=sha(fullrate_path),
                         raw_input_path=str(raw_path), raw_input_sha256=sha(raw_path),
                         output=artifact, mathematical_operator="native MF -> complex demod -> zero-phase circular FIR -> D decimation -> local interpolating DAS with (fc+offset) carrier",
                         reference_precision="FP64/complex128", additional_gain_fit=False,
                         baseline_unchanged=True)
            dump(dest / "comparison_to_v2.json", stats)
            audit.append(stats)
            print(f"{name}/{fixture}: relL2={stats.get('relative_l2')} rho={stats.get('magnitude_pearson')} pass={stats['passed']} seconds={stats['cpu_reference_seconds']:.2f}", flush=True)
    # Summarize all previously completed outputs, including resumed subsets.
    all_stats = []
    for name in CANDIDATES:
        for fixture in FIXTURE_NAMES:
            path = root / name / "cases" / fixture / "comparison_to_v2.json"
            if path.exists():
                all_stats.append(json.loads(path.read_text(encoding="utf-8")))
    manifest = dict(schema_version=3, shared_v2_bank=str(shared / "static"),
                    candidate_definitions=CANDIDATES, lut_phases=LUT_PHASES, kaiser_beta=KAISER_BETA,
                    numerical_static_assets=static_files, comparisons=all_stats,
                    full_suite_count=len(FIXTURE_NAMES)*len(CANDIDATES), completed_count=len(all_stats),
                    complete=len(all_stats)==len(FIXTURE_NAMES)*len(CANDIDATES),
                    exclusions=[dict(fixture="asymmetric_odd", reason="full-band template override fails subband precondition"),
                                dict(fixture="asymmetric_even", reason="full-band template override fails subband precondition")],
                    no_silent_fallback=True, no_threshold_relaxation=True)
    dump(root / "manifest.json", manifest)
    print(json.dumps(dict(output=str(root), completed=manifest["completed_count"], complete=manifest["complete"])), flush=True)


if __name__ == "__main__":
    main()
