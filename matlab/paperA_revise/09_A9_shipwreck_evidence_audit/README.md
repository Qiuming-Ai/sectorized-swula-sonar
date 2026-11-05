# A9 — Shipwreck-scene evidence audit

This historical source audit reads the three unchanged Exp-8 `.m.txt` snapshots in `legacy_sources/`. They are evidence inputs, not supported simulation entry points; no legacy MATLAB scripts, Field II, geometry data, or CUDA binaries are executed.

From the package root in MATLAB:

```matlab
addpath(fullfile(pwd,'matlab','paperA_revise','09_A9_shipwreck_evidence_audit'));
audit = run_A9_legacy_comparability_audit();
```

The default CSV is under `outputs/matlab/A9/`; pass a different output directory as the first argument to write elsewhere.

## Verdict

The STL geometry is useful as the input to a complex-scene stress test, but the existing outputs from the three legacy `PaperA_exp8_*_real_scene.m` scripts are **not a publishable baseline–CDM–FDM comparison**. The scripts do not hold the scene, aperture, waveform, noise realization, time-varying gain, or image reconstruction chain fixed. Consequently, the present images cannot isolate CDM waveform interference as the cause of image collapse.

The defensible description of a corrected experiment is **“STL-driven complex-scene channel simulation.”** The input is a triangular surface mesh processed by deterministic ray tracing, followed by a Born/single-scattering echo simulator. It is not measured sonar data and must not be described as a water-tank, transducer, shipwreck-survey, or sea experiment.

## Critical differences in the legacy scripts

| Control variable | Dense baseline | CDM | FDM | Consequence |
|---|---|---|---|---|
| Active scene | `model1 + model3` (ship + seabed) | `model2` only (`merged_model.stl`) | `model1 + model3` | CDM is not imaging the same object or background. |
| Secondary STL transform | scale 0.03, position `[0 15 8]` | scale 0.08, position `[10 15 8]` | scale 0.03, position `[0 15 8]` | Geometry, extent, bearing, and range distribution differ. |
| Receive array | 128 elements, about 0.5 lambda spacing | 32 elements, about 1 lambda spacing at 500 kHz | same as CDM | Does not match the revised manuscript's aperture-matched 32-element, 2-lambda SW-ULA. |
| Sparse transmit spacing | not applicable | lambda/4 | lambda/4 | Does not match the revised 12-element, 0.5-lambda, 18.3-mm transmit array. |
| Field of view | -65 to 65 deg | -60 to 60 deg | -60 to 60 deg | Does not match the revised common -50 to 50 deg FOV. |
| Pulse duration | 5 ms | about 5.03 ms (`N=503`, 100 kchip/s) | 0.5 ms | FDM time-bandwidth product is ten times smaller than declared in the manuscript. |
| Added noise | one measured-AWGN realization | four independently noisy copies summed | four independently noisy copies summed | The declared SNR and noise statistics are not common. |
| TVG | enabled, `K=10` | disabled | enabled, `K=40` | Range-dependent gain can create or hide apparent collapse. |
| Beamformer | fractional-delay (`das_plane_wave_fd`) | integer-delay (`das_plane_wave_id`) | fractional-delay (`das_plane_wave_fd`) | Reconstruction error is confounded with waveform error. |
| Decimation | 16 | 16 | 40 | Range grids and numerical bandwidth differ. |
| Manual image gain | none | none | one sector multiplied by 1.2 | Image intensity is not method-neutral. |
| Random control | ray seed only; seabed/noise not fixed before generation | same limitation | same limitation | A repeatable, paired comparison is unavailable. |

The source-level audit is reproduced by `run_A9_legacy_comparability_audit.m` and written to `outputs/matlab/A9/A9_legacy_comparability_audit.csv` by default.

## Minimum protocol before the result can enter the paper

1. Generate the ray-traced scatterer list once and save `point_position`, `point_amplitudes`, mesh identity, transform, and random seed. Load this identical list for all three methods.
2. Use the manuscript geometry: dense 128-element receive ULA at 0.5 lambda; proposed 32-element receive SW-ULA at 2 lambda; both apertures approximately 0.21 m; proposed 12-element transmit ULA at 0.5 lambda.
3. Use the common -50 to 50 deg FOV and the six revised steering angles `[-41.25, -24.25, -8, 8, 24.25, 41.25]` deg.
4. Use 5-ms pulses in all methods. Use the declared 499-chip ZC configuration for CDM and six 20-kHz LFM subbands within the 100-kHz allocation for FDM.
5. Normalize transmit energy or state an explicit equal-input-SNR image-quality comparison. Add exactly one paired noise realization per method and seed; do not sum four noisy copies.
6. Use the same receive-aperture taper, TVG rule, fractional-delay beamformer, output range grid, angular grid, display normalization, and dynamic range.
7. Export complex images before plotting. Derive any foreground/support mask from the dense reference once, and apply the identical mask to all methods.
8. Report amplitude-aligned NMSE, global SSIM, support leakage, and separate angular/range spreading. Repeat the noisy comparison for at least 30 paired seeds and report confidence intervals.

## Appropriate evidential strength

- Current legacy figures: exploratory debugging evidence only.
- Corrected common-scene channel-level rerun: strong numerical evidence that finite-pulse CDM interference accumulates in a geometrically complex, multi-range scene.
- Neither version: real-acoustic validation.

The corrected A11 channel-level rerun is now complete. It retains the dense-baseline script's ship scale, pose, placement, coordinate permutation, 3-D seabed, and one-bounce CUDA-visible arrival cloud. Because the receive array is one-dimensional along x, the image angle is the direction-cosine coordinate `asin(x/r)`, with `r=sqrt(x^2+y^2+z^2)`; using `atan2(x,z)` would incorrectly discard unresolved elevation. The old 45-deg direction is preview metadata and does not define the image axis when directional gains are disabled. One common element-level channel and reconstruction chain is then applied to dense, CDM, and FDM. The remaining angle--slant-range curvature reflects elevation collapse and is not removed by ground-truth-aided warping. The result is strong numerical evidence of method-dependent behavior in a geometrically complex scene, but it is not real-acoustic, elevation-resolved, multipath, or reverberant validation.
