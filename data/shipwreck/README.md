# Shipwreck data

The author generated the STL and its derived geometry. The STL, simulated arrival catalogue, and geometry-derived image arrays are released under [Apache-2.0](../../LICENSE). They are numerical-model data, not measured underwater echoes or an as-surveyed wreck.

| File in the companion asset | Role |
|---|---|
| `ship.stl` | Original simulated-scene mesh. |
| `A11_raytrace_arrivals.csv` | Archived arrival catalogue including coordinates, range, azimuth, elevation, amplitude and phase; input to the numerical rerun. |
| `A11_plot_data.mat` | Lossless export of the three reference arrays used by Figure 13: `range_m`, `angle_deg`, and `truth_magnitude`. Unused MATLAB workspace objects are omitted. |
| `noisy_common_per_element_noise_10dB.npz` | Original saved complex reconstruction grids (`dense`, `cdm`, `fdm`) and coordinates (`ranges_m`, `angles_deg`) for the representative common-noise 10 dB example. |

Download `shipwreck-data_v10.zip` from the release that contains this repository version. Extract the archive at the repository root; its paths begin with `data/shipwreck/geometry/`. Do not commit the extracted directory to Git. From the repository root:

```sh
python scripts/check_release.py --shipwreck
python code/scenes/plot_shipwreck.py
```

The final Python rerun in `code/scenes/run_final_scenes.py --suite ship` starts from the archived arrival catalogue. It does not rebuild it from STL. See the [MATLAB guide](../../matlab/README.md) for the historical geometry-generation chain. The `all_trials.csv` and `summary_mean_ci.csv` files tracked in Git describe 30 phase/noise realizations of the same fixed geometry. They support Figure 14 and Table 7 independently of the companion asset.
