"""Compute reproducible image-quality metrics for the reviewer-ready figures.

The script evaluates the six ``main`` data sets at input SNR = -15 dB.  All
methods use the same ground-truth masks and search windows.  Metrics are
computed on linear image amplitude before the display dynamic-range clip.

Definitions
-----------
* localization error: mean absolute offset of the strongest sample inside a
  +/-2 deg, +/-0.04 m window around each ground-truth target;
* CNR: (mean target amplitude - mean background amplitude) divided by the
  background standard deviation;
* support leakage: fraction of image energy outside the union of target guard
  regions (+/-2.5 deg, +/-0.08 m);
* maximum false response: strongest sample outside all guard regions relative
  to the global image peak;
* false-peak count: local maxima outside all guard regions above -10 dB,
  using a 1 deg x 0.015 m non-maximum-suppression neighbourhood.

The dense target cluster makes these deliberately conservative scene-level
metrics; point-spread width and PSLR/ISLR are reported separately by A2.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.ndimage import maximum_filter


METHODS = ("dense", "cdm", "fdm", "random", "nested", "coprime")
METHOD_LABELS = {
    "dense": "Dense ULA",
    "cdm": "CDM SW-ULA",
    "fdm": "FDM SW-ULA",
    "random": "Random sparse",
    "nested": "Nested sparse",
    "coprime": "Coprime sparse",
}

LOCALIZATION_HALF_ANGLE_DEG = 2.0
LOCALIZATION_HALF_RANGE_M = 0.04
CORE_HALF_ANGLE_DEG = 0.6
CORE_HALF_RANGE_M = 0.01
GUARD_HALF_ANGLE_DEG = 2.5
GUARD_HALF_RANGE_M = 0.08
FALSE_PEAK_THRESHOLD_DB = -10.0
FALSE_PEAK_NMS_ANGLE_DEG = 1.0
FALSE_PEAK_NMS_RANGE_M = 0.015


def _load(path: Path) -> tuple[np.ndarray, ...]:
    raw = loadmat(path, squeeze_me=True)
    angle = np.asarray(raw["angle_deg"], dtype=float).reshape(-1)
    ranges = np.asarray(raw["range_m"], dtype=float).reshape(-1)
    image = np.abs(np.asarray(raw["image_linear"], dtype=float))
    if image.shape == (angle.size, ranges.size):
        image = image.T
    if image.shape != (ranges.size, angle.size):
        raise ValueError(f"Unexpected image shape in {path}: {image.shape}")
    image /= np.nanmax(image)
    target_angle = np.asarray(raw["target_angle_deg"], dtype=float).reshape(-1)
    target_range = np.asarray(raw["target_range_m"], dtype=float).reshape(-1)
    return angle, ranges, image, target_angle, target_range


def _union_mask(
    angle: np.ndarray,
    ranges: np.ndarray,
    target_angle: np.ndarray,
    target_range: np.ndarray,
    half_angle: float,
    half_range: float,
) -> np.ndarray:
    mask = np.zeros((ranges.size, angle.size), dtype=bool)
    for theta, radius in zip(target_angle, target_range):
        mask |= (
            (np.abs(angle[None, :] - theta) <= half_angle)
            & (np.abs(ranges[:, None] - radius) <= half_range)
        )
    return mask


def _localization_errors(
    angle: np.ndarray,
    ranges: np.ndarray,
    image: np.ndarray,
    target_angle: np.ndarray,
    target_range: np.ndarray,
) -> tuple[float, float, float, float]:
    angle_errors: list[float] = []
    range_errors: list[float] = []
    for theta, radius in zip(target_angle, target_range):
        angle_idx = np.flatnonzero(
            np.abs(angle - theta) <= LOCALIZATION_HALF_ANGLE_DEG
        )
        range_idx = np.flatnonzero(
            np.abs(ranges - radius) <= LOCALIZATION_HALF_RANGE_M
        )
        patch = image[np.ix_(range_idx, angle_idx)]
        local_r, local_a = np.unravel_index(np.nanargmax(patch), patch.shape)
        angle_errors.append(abs(float(angle[angle_idx[local_a]]) - theta))
        range_errors.append(abs(float(ranges[range_idx[local_r]]) - radius))
    return (
        float(np.mean(angle_errors)),
        float(np.max(angle_errors)),
        float(np.mean(range_errors)),
        float(np.max(range_errors)),
    )


def _false_peak_count(
    angle: np.ndarray,
    ranges: np.ndarray,
    image: np.ndarray,
    background_mask: np.ndarray,
) -> int:
    angle_step = float(np.median(np.diff(angle)))
    range_step = float(np.median(np.diff(ranges)))
    nms_angle = max(3, int(round(FALSE_PEAK_NMS_ANGLE_DEG / angle_step)))
    nms_range = max(3, int(round(FALSE_PEAK_NMS_RANGE_M / range_step)))
    local_max = maximum_filter(image, size=(nms_range, nms_angle), mode="nearest")
    threshold = 10.0 ** (FALSE_PEAK_THRESHOLD_DB / 20.0)
    peaks = background_mask & (image >= threshold) & np.isclose(image, local_max)
    return int(np.count_nonzero(peaks))


def compute_metrics(data_dir: Path) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for method in METHODS:
        path = data_dir / f"main_{method}_snr_m15.mat"
        angle, ranges, image, target_angle, target_range = _load(path)

        # Restrict evaluation to the displayed reviewer figure.
        crop = (
            (angle >= -50.0),
            (angle <= 50.0),
            (ranges >= 4.80),
            (ranges <= 5.70),
        )
        angle_idx = np.flatnonzero(crop[0] & crop[1])
        range_idx = np.flatnonzero(crop[2] & crop[3])
        angle_eval = angle[angle_idx]
        range_eval = ranges[range_idx]
        image_eval = image[np.ix_(range_idx, angle_idx)]

        core = _union_mask(
            angle_eval,
            range_eval,
            target_angle,
            target_range,
            CORE_HALF_ANGLE_DEG,
            CORE_HALF_RANGE_M,
        )
        guard = _union_mask(
            angle_eval,
            range_eval,
            target_angle,
            target_range,
            GUARD_HALF_ANGLE_DEG,
            GUARD_HALF_RANGE_M,
        )
        background = ~guard

        target_mean = float(np.mean(image_eval[core]))
        background_mean = float(np.mean(image_eval[background]))
        background_std = float(np.std(image_eval[background], ddof=1))
        cnr = (target_mean - background_mean) / background_std
        total_energy = float(np.sum(np.square(image_eval)))
        leakage = float(np.sum(np.square(image_eval[background])) / total_energy)
        false_level_db = 20.0 * np.log10(
            float(np.max(image_eval[background])) / float(np.max(image_eval))
        )
        mean_a, max_a, mean_r, max_r = _localization_errors(
            angle_eval,
            range_eval,
            image_eval,
            target_angle,
            target_range,
        )

        rows.append(
            {
                "method": method,
                "method_label": METHOD_LABELS[method],
                "mean_angle_error_deg": mean_a,
                "max_angle_error_deg": max_a,
                "mean_range_error_m": mean_r,
                "max_range_error_m": max_r,
                "cnr_linear": cnr,
                "support_leakage": leakage,
                "max_false_response_db": false_level_db,
                "false_peak_count_above_m10db": _false_peak_count(
                    angle_eval, range_eval, image_eval, background
                ),
            }
        )
    return rows


def main() -> None:
    root = Path(__file__).resolve().parent
    data_dir = root / "data"
    output = root / "results" / "A8_quantitative_metrics.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = compute_metrics(data_dir)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row)
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
