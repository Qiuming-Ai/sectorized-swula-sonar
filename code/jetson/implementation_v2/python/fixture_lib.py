"""Binary contract and independently implemented FP64 discrete reference.

The reference accepts quantized input, but not a fitted gain, GPU output, or a
reconstruction from the frozen model. Only the emitter/template construction is
shared with that frozen model. All arrays use little-endian, C-major storage.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy import signal

MODEL_SHA256 = "4b280863495d4d0f677f36dac79736cedd9d2ba009a3cacf41849cfce7db8f64"
METHODS = ("dense", "cdm", "fdm")
FIXED_GATES = dict(stage_relative_l2_max=1e-4, image_relative_l2_max=1e-3,
                   magnitude_pearson_min=.999, zero_max_abs=1e-7,
                   peak_range_difference_max_m=.0025,
                   peak_angle_difference_max_deg=.25)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def dump(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
        allow_nan=False, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + "\n", encoding="utf-8")

def read_meta(path):
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        key, val = line.split("=", 1)
        try:
            out[key] = float(val) if any(c in val for c in ".eE") else int(val)
        except ValueError:
            out[key] = val
    return out

def read_method(path):
    path = Path(path)
    meta = read_meta(path / "meta.txt")
    m, n, k = int(meta["channels"]), int(meta["fft_samples"]), int(meta["waveforms"])
    return dict(meta=meta,
        g=np.fromfile(path / "G64.bin", dtype="<c16").reshape(n),
        q=np.fromfile(path / "Q64.bin", dtype="<c16").reshape(k, n),
        positions=np.fromfile(path / "positions.bin", dtype="<f8").reshape(m),
        angles=np.fromfile(path / "angles.bin", dtype="<f8"),
        ranges=np.fromfile(path / "ranges.bin", dtype="<f8"),
        sector_map=np.fromfile(path / "sector_map.bin", dtype="<i4"))

def read_raw(path, meta):
    return np.fromfile(path, dtype="<c8").reshape(int(meta["channels"]), int(meta["record_samples"]))

def receive_fft(raw, bank):
    """FE spectral stage: one common receive filter, exactly once."""
    n = int(bank["meta"]["fft_samples"])
    return np.fft.fft(np.asarray(raw, np.complex128), n=n, axis=1) * bank["g"][None, :]

def matched_time(fe, bank, waveform):
    """Unshifted discrete correlation: ifft(FE * conjugate(q))."""
    return np.fft.ifft(fe * bank["q"][waveform].conj()[None, :], axis=1)

def independent_image(raw, bank, return_stages=False):
    meta = bank["meta"]
    n, os = int(meta["fft_samples"]), int(meta["oversampling"])
    m, kmax = int(meta["channels"]), int(meta["waveforms"])
    r, a = bank["ranges"], bank["angles"]
    wr = np.hamming(m); wr /= wr.sum()
    # Exact equal-sector rule, including final endpoint, same as the paper model.
    sector = bank["sector_map"]
    fe = receive_fft(raw, bank)
    image = np.zeros((len(r), len(a)), dtype=np.complex128)
    stages = [] if return_stages else None
    for k in range(kmax):
        ix = np.arange(len(a)) if kmax == 1 else np.flatnonzero(sector == k)
        mf = matched_time(fe, bank, k)
        if return_stages:
            stages.append(mf)
        for channel, position in enumerate(bank["positions"]):
            # scipy complex Fourier resampling splits even-length Nyquist into
            # positive/negative halves. This convention is part of the contract.
            fine = signal.resample(mf[channel], n * os)
            delay = (2 * (r[:, None] - float(meta["gate_range"])) - position * np.sin(a[ix])[None, :]) / float(meta["c"])
            index = delay * float(meta["fs"]) * os
            low = np.floor(index).astype(np.int64)
            if np.any(low < 0) or np.any(low + 1 >= n * os):
                raise ValueError("Declared image grid lies outside guarded linear record")
            fraction = index - low
            samples = fine[low] * (1 - fraction) + fine[low + 1] * fraction
            image[:, ix] += wr[channel] * samples * np.exp(2j * np.pi * float(meta["fc"]) * delay)
    return (image, fe, stages) if return_stages else image

def comparison(reference, candidate, bank=None, kind="image"):
    ref, cand = np.asarray(reference), np.asarray(candidate)
    if ref.shape != cand.shape:
        raise ValueError(f"Shape mismatch {ref.shape} != {cand.shape}")
    diff = cand - ref
    scale = np.linalg.norm(ref.ravel())
    out = dict(shape=list(ref.shape), max_abs_error=float(np.max(np.abs(diff))),
               reference_l2=float(scale), candidate_l2=float(np.linalg.norm(cand.ravel())),
               fitted_gain_used=False, all_finite=bool(np.isfinite(cand).all()))
    if scale == 0:
        out.update(relative_l2=None, zero_max_abs=float(np.max(np.abs(cand))),
                   passed=bool(np.max(np.abs(cand)) <= FIXED_GATES["zero_max_abs"] and out["all_finite"]))
        return out
    relative = float(np.linalg.norm(diff.ravel()) / scale)
    out["relative_l2"] = relative
    limit = FIXED_GATES["image_relative_l2_max"] if kind == "image" else FIXED_GATES["stage_relative_l2_max"]
    passed = relative <= limit and out["all_finite"]
    if kind == "image":
        u, v = np.abs(ref).ravel(), np.abs(cand).ravel()
        u -= u.mean(); v -= v.mean()
        denominator = np.linalg.norm(u) * np.linalg.norm(v)
        rho = float(np.dot(u, v) / denominator) if denominator else None
        out["magnitude_pearson"] = rho
        if rho is not None:
            passed &= rho >= FIXED_GATES["magnitude_pearson_min"]
        elif np.linalg.norm(u) > 0:
            # Only a constant reference makes correlation not applicable.
            # A constant candidate cannot pass against a varying reference.
            passed = False
        if bank is not None:
            p0 = np.unravel_index(np.argmax(np.abs(ref)), ref.shape)
            p1 = np.unravel_index(np.argmax(np.abs(cand)), cand.shape)
            rd = abs(float(bank["ranges"][p1[0]] - bank["ranges"][p0[0]]))
            ad = abs(float(np.rad2deg(bank["angles"][p1[1]] - bank["angles"][p0[1]])))
            out.update(peak_range_difference_m=rd, peak_angle_difference_deg=ad,
                       reference_peak_index=list(map(int, p0)), candidate_peak_index=list(map(int, p1)))
            passed &= rd <= FIXED_GATES["peak_range_difference_max_m"] + 1e-12
            passed &= ad <= FIXED_GATES["peak_angle_difference_max_deg"] + 1e-12
    out["passed"] = bool(passed)
    return out
