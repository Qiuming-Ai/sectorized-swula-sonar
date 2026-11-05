# Manuscript-to-artifact map

Numbering follows the v08 manuscript that supplies the archived numerical evidence for [the article](https://doi.org/10.1016/j.rineng.2026.113776), not historical filenames. “Replot” means plotting the stored numerical inputs, not necessarily byte-identical PDF output. An explanatory schematic has no experimental dataset.

| Figure | Evidence in this release | Boundary |
|---|---|---|
| 1–2 | [Design notes and source mapping](../code/design/README.md) | Explanatory aliasing/masking graphics; distinguish exact historical art from a model illustration. |
| 3 | `data/spatial/Figure03_transmit_beampattern_data.csv`; [design scripts](../code/design/README.md) | Six 450-kHz curves, not the full broadband search. |
| 4(a,b) | [Legacy design-waveform examples](../code/design/README.md) | Ideal/zero-Doppler design examples, not the final front-end-aware Figure 5. |
| 5 | `data/waveforms/figure05/`; `code/plots/figure05/plot_figure5_legend.py`; `data/waveforms/ambiguity/` | Saved full/local plotted pairs, spectra, and worst/ordered-pair tables; large complete correlation bank is omitted. |
| 6–7 | [Design and pipeline notes](../code/design/README.md); `code/common/` | Architecture schematics; no field data. |
| 8 | `data/spatial/A14_angular_profiles.csv`; [design scripts](../code/design/README.md) | Receiver-only and assigned-waveform angular products. |
| 9 | `data/waveforms/two_way/`; `code/wave_spatial/plot_final_figures.py`; `data/spatial/` | Full-bank vs centre-frequency spatial approximation, not a dense-scene reconstruction. |
| 10 | `data/point30m/point30m_main_minus15dB_arrays.npz`; `code/scenes/plot_submission_figures.py` | Common per-element noise, −15 dB Dense-reference input SNR, one paired realization. |
| 11 | `data/point30m/`, `data/dense1000x50/`; `code/scenes/plot_submission_figures.py` | Thirty noise trials per method, SNR and noise rule for each fixed scene. |
| 12 | `data/density/`; `code/plots/plot_density_main.py`; `code/density_snapshot/` | 99 displayed rows; 11 conditions × three generated scenes × three methods. |
| 13 | [Shipwreck release data](../data/shipwreck/README.md); `code/scenes/plot_shipwreck.py`; `code/scenes/run_final_scenes.py --suite ship` | Replot from image arrays or rerun from the frozen arrival catalogue. The [historical MATLAB/CUDA chain](../matlab/paperA_revise/11_A11_shipwreck_common_scene/README.md) separately supports STL-to-arrivals generation with a compatible local MEX toolchain. |
| 14 | `data/shipwreck/`; `code/plots/redraw_figure14_bars.py` | Scalar metrics only; phase/noise variation of one geometry. |
| 15 | `data/jetson/`; [Jetson entry points](../code/jetson/README.md) | Archived benchmark logs; see the Jetson notes for the distinction between summary bars and the raw timing distribution. |

| Table or section | Evidence |
|---|---|
| Table 1: simulation parameters | `code/common/numerical_model.py` (`Config` and waveform constants), `data/waveforms/figure05/configuration.json`, `code/scenes/run_final_scenes.py`, and `data/inputs/`. |
| Section 2 coverage and grating constraints | `data/coverage/` contains the candidate search, implemented six-sector audit and grid convergence. `code/coverage/run_coverage.py` calls the frozen coverage implementation. The accepted result is limited to the tested candidate family and thresholds. |
| Table 2: centre-frequency spatial products | [A14 design notes](../code/design/README.md) and `data/spatial/A14_angular_profiles.csv`. This is not the Figure 9 all-bank response. |
| Table 3: envelope controls | `data/waveforms/window/window_template_metrics.csv`; `code/wave_spatial/run_wave_spatial.py --only window`. |
| Table 4: point detection/localization | `data/point30m/all_trials.csv`, `summary_mean_ci.csv`, independent noise calibration and false-alarm files. Select common per-element noise at −15 dB. |
| Table 5: NMS sensitivity | `data/nms_sensitivity/`; `code/scenes/nms_sensitivity.py` with the full independent test-noise cache. The original threshold is held fixed. |
| Table 6: DAS/OMP | `data/das_omp/` and `code/common/competitive.py --min-spacing-width`; random-array comparison within paired measurements. |
| Table 7: shipwreck results | `data/shipwreck/summary_mean_ci.csv` and `all_trials.csv`, common per-element noise at 10 dB. |
| Tables 8–9: platform and runtime | [Jetson build and measurement notes](../code/jetson/README.md) and `data/jetson/`. Platform descriptors are not new measurements. |
| Table 10: trade-offs | References the envelope, finite-pulse ambiguity, shipwreck and Jetson datasets above; no independent run. |

The public bundle is based on frozen evidence, not the article-submission folder. The publication-stage MATLAB scripts, final Python simulation snapshots, and Jetson replay implementation are not interchangeable.

The [desktop MATLAB timing archive](../data/desktop-timing/README.md) is separate from the measured Jetson workload in Table 9 and Figure 15.
