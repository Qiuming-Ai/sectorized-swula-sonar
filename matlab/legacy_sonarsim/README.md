# Historical SonarSim experiments and desktop MATLAB timing

This module preserves the 28 paper-specific `SonarSim/PaperA_*.m` entries, the directly required author-side MATLAB/CUDA helpers, and their small figure and timing records. The historical scripts are not substituted for the final manuscript's uniform numerical suite or the separate Jetson experiment.

## Verify the recorded desktop result

From the repository root:

```matlab
addpath(fullfile(pwd,'matlab','legacy_sonarsim'));
legacy_sonarsim_smoke();
```

This needs no external HDF5, STL, or CUDA MEX. The 100-row archived `data/TimingData_20251126_114034.csv` stores seconds, and gives Dense / CDM / FDM means of **265.608743 / 60.105308 / 31.266306 ms**, matching **265.61 / 60.11 / 31.27 ms** in the manuscript's desktop MATLAB paragraph. Do not substitute the later timing CSV or the Jetson CSVs.

## Optional exact-input replay

Extract the optional `desktop-timing-data_v10.zip` release asset so the three HDF5 files lie directly under `data/desktop-timing/overlay/`. [Input requirements](data/input_requirements.json) gives filenames, sizes, SHA-256 values, and dimensions. The files total 110,161,193 bytes and are not in the source tree.

| File | Samples/channel | Channels | Pings | Stored SNR / velocity |
|---|---:|---:|---:|---|
| `BaselineData_20251121_150913.h5` | 712 | 128 | 100 | −20 dB / 0.2 m/s |
| `CDMData_20251125_160029.h5` | 717 | 32 | 100 | +10 dB / 0 m/s |
| `FDMData_20251126_112556.h5` | 713 | 32 | 100 | −10 dB / 0.2 m/s |

These are generated complex-baseband data, not acoustic acquisitions. The files have `/raw_data/ping_1` through `ping_100` groups with `real` and `imag` float32 datasets, and `/raw_data/.attributes` for geometry, sampling and matching templates. MATLAB reads each ping as samples × channels. The FDM metadata also stores its 85-tap second filter and decimation factor 3. The local data-maker sources document the synthetic generation lineage, but have changed since the archived inputs were created; they are not a byte-identical regeneration recipe.

```matlab
addpath(fullfile(pwd,'matlab','legacy_sonarsim'));
input_dir = fullfile(pwd,'data','desktop-timing','overlay');
legacy_sonarsim_smoke(input_dir);          % one CPU ping per method
new_csv = [tempname(tempdir), '.csv'];     % never replace the archived CSV
T = CollectTimingData(input_dir,new_csv); % full 100-ping replay per method
```

The collector preserves the historical Dense → FDM → CDM order and timed matching/DAS blocks. FDM includes its subband mix, second FIR and resampling inside the timer; HDF5 reads and the optional historical post-timer FDM plot are outside. A third argument `true` enables that plot. Replay numbers depend on MATLAB, hardware and process state; the recorded table is reproduced from the CSV, not from a new stopwatch. These desktop workloads differ in raw input sizes and conditions and are not a matched desktop/Jetson speedup experiment.

## Historical experiment map

| Entry family | Role |
|---|---|
| `PaperA_exp0`, `exp1*`, `exp2` | Early point-scatterer, receive-processing and subband prototypes |
| `PaperA_exp3_*` | Dense/CDM/FDM noise and motion prototypes, including the original baseband entries used by the A8 exporter |
| `PaperA_exp4_*` | Random/nested/coprime legacy controls; also read by A8 |
| `PaperA_exp5_*`, `plot_paperA_exp5_*`, `test_plot_paperA_exp5*` | Transmit-pattern/pulse-compression examples; three Exp5 inputs and one Exp7 input are stored under `data/figures/` |
| `PaperA_exp6_*`, `PaperA_exp7_0`, `plot_paperA_exp7` | Alternative imaging cases and archived angular HPBW/PSLR analysis |
| `PaperA_exp8_*` | Historical mesh-scatterer prototypes; not the final common-scene comparison |
| `CollectTimingData`, `*ProcessFD`, `*DataMaker` | Desktop timing replay and synthetic fixture-generation lineage |

Use MATLAB with Signal Processing Toolbox for FIR/chirp/resampling and Communications Toolbox for experiments using `awgn`. Run legacy generators in a scratch working directory after adding this module to the MATLAB path; some generate figures or data in the working directory. The archived figure MAT files contain only numerical plot inputs. The package intentionally excludes unrelated projects, MATLAB toolbox files, autosaves and GUI utilities.

### Optional CUDA and mesh paths

The CDM baseband Exp3/A8 chain uses the optional echo MEX. The Exp8 scripts additionally need the CUDA ray tracer. With a MATLAB-compatible NVIDIA toolkit and host compiler, use `build_legacy_mex("echo")` or `build_legacy_mex("all")`. The source is included; no MEX binaries or SDK files are distributed. Build outputs go to ignored `build_mex/`. The packaging host's MATLAB-bundled CUDA headers rejected its installed MSVC (C1189) in the ray-tracer build check. No compiler-check bypass was used, and no live CUDA execution is claimed for this legacy module.

Exp8 resolves the exact `ship.stl` first from extracted `data/shipwreck/geometry/`, then from the module's local `Models/` directory. The exact `merged_model.stl` is a separate opt-in local asset and is not redistributed. See the input manifest for both hashes; do not substitute a different mesh. All three original scripts read both meshes, although Dense/FDM use ship plus seabed and the CDM prototype uses the merged object. That difference is preserved and is why Exp8 is not a fair method comparison.

## Source provenance and packaging checks

Original source locations are `SonarSim/` (entry scripts, timing CSV, and figure MATs) and `SonarSim/Lib/` (direct helpers). The release copy limits path setup to this module, removes unused broad workspace and Field II path additions, changes the three stale FD-DAS callers to the local two-output interface, and makes timing, mesh, and figure paths explicit. Timing loops, signal and geometry parameters, and A8 numerical text hooks are preserved. The replay collector avoids overwriting outputs; DataMaker and Exp7 saves reject an existing destination.

Native MATLAB R2025b smoke passed: archived means, local FD-DAS shape/finite output, four readable figure MATs, and one CPU ping per method using the exact original HDF5 files read-only. The full archived GPU/STL experiments were not run during release preparation.
