"""Independent FP64 reduced-rate matched-filter / local-delay reference.

This module never reads CUDA outputs. The input bank and complex-FP32 raw
fixtures are frozen v2 assets. Downsampling candidates do NOT modify that bank.
The only common implementation detail is the preregistered mathematical LUT
definition, explicitly exported in FP64 and quantized FP32 variants.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy import signal

HERE = Path(__file__).resolve().parent
V2_PYTHON = HERE.parents[1] / "implementation_v2" / "python"
sys.path.insert(0, str(V2_PYTHON))
from fixture_lib import read_method, read_raw, receive_fft, matched_time, comparison, sha, dump

CANDIDATES = {
    "d4_linear": dict(decimation=4, fir_taps=129, cutoff_hz=30000., interpolation_taps=0),
    "d4_sinc16": dict(decimation=4, fir_taps=129, cutoff_hz=30000., interpolation_taps=16),
    "d8_linear": dict(decimation=8, fir_taps=257, cutoff_hz=17500., interpolation_taps=0),
    "d8_sinc32": dict(decimation=8, fir_taps=257, cutoff_hz=17500., interpolation_taps=32),
}
OFFSETS = np.array([-24000., -8000., 40000., 24000., 8000., -38000.], dtype=np.float64)
LUT_PHASES = 1024
KAISER_BETA = 8.6
FIXTURE_NAMES = ["paper_clean"] + [f"noise{k:02d}" for k in range(20)] + [
    "zero", "boundary_points", "offgrid", "complex_impulse"]


def ensure_applicable(bank):
    """Fail fast rather than silently applying a different full-band model."""
    meta = bank["meta"]
    expected = {"channels": 32, "record_samples": 8192, "fft_samples": 16384,
                "waveforms": 6, "fs": 400000., "fc": 450000., "c": 1500.}
    for key, value in expected.items():
        if meta.get(key) != value:
            raise ValueError(f"Unsupported frozen FDM bank: {key}={meta.get(key)}, expected {value}")
    if meta.get("unit_template_override", 0):
        raise ValueError("Full-band asymmetric-template override is not a subband FDM input")
    if len(bank["sector_map"]) != len(bank["angles"]):
        raise ValueError("Sector-map length differs from angular grid")


def interpolation_lut(taps, phases=LUT_PHASES):
    if taps == 0:
        return np.empty((phases + 1, 0), dtype=np.float64)
    shifts = np.arange(-taps // 2 + 1, taps // 2 + 1)
    distance = np.arange(phases + 1)[:, None] / phases - shifts[None, :]
    radius = taps / 2
    inside = np.abs(distance) < radius
    window = np.zeros_like(distance)
    window[inside] = np.i0(KAISER_BETA * np.sqrt(1 - (distance[inside] / radius) ** 2)) / np.i0(KAISER_BETA)
    weights = np.sinc(distance) * window
    weights /= weights.sum(axis=1, keepdims=True)
    return weights


def candidate_bank(bank, name):
    ensure_applicable(bank)
    cfg = CANDIDATES[name]
    meta = bank["meta"]
    n, fs = int(meta["fft_samples"]), float(meta["fs"])
    taps = signal.firwin(cfg["fir_taps"], cfg["cutoff_hz"], fs=fs,
                        window=("kaiser", KAISER_BETA), scale=True)
    delay = (cfg["fir_taps"] - 1) / 2
    frequency = np.fft.fftfreq(n, 1 / fs)
    h = np.fft.fft(taps, n) * np.exp(2j * np.pi * frequency * delay / fs)
    mixing = np.exp(-2j * np.pi * OFFSETS[:, None] * np.arange(n)[None, :] / fs)
    return dict(name=name, **cfg, h=h, mix=mixing, offsets=OFFSETS.copy(),
                lut=interpolation_lut(cfg["interpolation_taps"]), lut_phases=LUT_PHASES,
                fir_delay_compensated_samples=delay, fs_low=fs / cfg["decimation"],
                low_samples=n // cfg["decimation"])


def filtered_decimated(matched, candidate, k):
    mixed = matched * candidate["mix"][k][None, :]
    return np.fft.ifft(np.fft.fft(mixed, axis=1) * candidate["h"][None, :], axis=1)[:, ::candidate["decimation"]]


def local_interpolate(record, indices, candidate):
    low = np.floor(indices).astype(np.int64)
    frac = indices - low
    taps = candidate["interpolation_taps"]
    if not taps:
        if low.min() < 0 or low.max() + 1 >= len(record):
            raise ValueError("Linear delayed sample lies outside guarded record")
        return record[low] * (1 - frac) + record[low + 1] * frac
    shifts = np.arange(-taps // 2 + 1, taps // 2 + 1)
    if low.min() + shifts[0] < 0 or low.max() + shifts[-1] >= len(record):
        raise ValueError("Sinc delayed sample lies outside guarded record")
    phase = frac * candidate["lut_phases"]
    phase_lo = np.floor(phase).astype(np.int64)
    phase_blend = phase - phase_lo
    lut = candidate["lut"]
    samples = np.zeros(indices.shape, dtype=np.complex128)
    # Loop over tap rather than constructing an enormous [R,A,tap] tensor.
    for j, shift in enumerate(shifts):
        weight = lut[phase_lo, j] * (1 - phase_blend) + lut[phase_lo + 1, j] * phase_blend
        samples += record[low + shift] * weight
    return samples


def image_from_mf(matched_records, bank, candidate, return_low=False):
    """DAS on genuine low-rate records; no full-record Fourier upsampling."""
    meta = bank["meta"]
    r, a = bank["ranges"], bank["angles"]
    wr = np.hamming(int(meta["channels"])); wr /= wr.sum()
    image = np.zeros((len(r), len(a)), dtype=np.complex128)
    low_records = [] if return_low else None
    for k, mf in enumerate(matched_records):
        ix = np.flatnonzero(bank["sector_map"] == k)
        low = filtered_decimated(mf, candidate, k)
        if return_low:
            low_records.append(low)
        for channel, position in enumerate(bank["positions"]):
            delay = (2 * (r[:, None] - float(meta["gate_range"])) - position * np.sin(a[ix])[None, :]) / float(meta["c"])
            samples = local_interpolate(low[channel], delay * candidate["fs_low"], candidate)
            image[:, ix] += wr[channel] * samples * np.exp(2j * np.pi * (float(meta["fc"]) + candidate["offsets"][k]) * delay)
    return (image, low_records) if return_low else image


def independent_image(raw, bank, candidate, return_stages=False):
    ensure_applicable(bank)
    fe = receive_fft(raw, bank)
    matched = [matched_time(fe, bank, k) for k in range(int(bank["meta"]["waveforms"]))]
    result = image_from_mf(matched, bank, candidate, return_low=return_stages)
    if return_stages:
        image, low = result
        return image, fe, matched, low
    return result


def save_binary(path, array, dtype):
    """Binary numerical artifacts, not source-file editing."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.ascontiguousarray(array, dtype=np.dtype(dtype))
    arr.tofile(path)
    return dict(path=str(path), dtype=arr.dtype.str, shape=list(arr.shape),
                bytes=path.stat().st_size, sha256=sha(path))
