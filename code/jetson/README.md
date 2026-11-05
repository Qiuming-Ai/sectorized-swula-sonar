# Jetson receive-reconstruction benchmark

This directory supports two different tasks: reanalyzing the **recorded paper results** on an ordinary computer, and building a **new replay** on Jetson. The historical source and data are not overwritten by either task.

## Recalculate the paper results

Run from the repository root with Python 3.11 or later. The verifier uses only the Python standard library and does not need a GPU, CUDA, MATLAB, or network access.

```sh
python code/jetson/verify_records.py
```

The verifier checks seven frozen source files and the frozen CPU model hash. It recalculates all 1,800 measured wall times, nine explicit device-memory ledgers, and nine energy integrals directly from timestamped current VDD_IN readings. It checks the 85,439 completed blocks against the power-campaign block logs and verifies the historical 75-case own-reference receipt. It does not claim to rerun those historical GPU correctness checks.

| Configuration | Time (ms), mean ± sample SD | P95 (ms) | Explicit CUDA high-water (MiB) | Module power (W), mean ± sample SD | Module energy (J/block), mean ± sample SD |
|---|---:|---:|---:|---:|---:|
| Dense LFM | 20.476 ± 0.096 | 20.613 | 282.022 | 17.314 ± 0.053 | 0.35592 ± 0.00123 |
| CDM-ZC | 16.431 ± 0.085 | 16.549 | 72.059 | 19.369 ± 0.032 | 0.31911 ± 0.00042 |
| FDM-LFM | 7.773 ± 0.075 | 7.911 | 14.009 | 16.396 ± 0.005 | 0.12812 ± 0.00010 |

These match the manuscript table labeled `tab:runtime` (Table 9 in the standalone v08 package). Each method contributes 600 complete timed calls across three sessions. The unit for power/energy variation is an independent workload window, with three windows per method. P95 uses type-7 interpolation. The memory figures include explicitly allocated shared cuFFT workspace; system RAM and process RSS are separate observations and are not added to device allocation.

Optional figures require the analysis dependencies:

```sh
python -m pip install -r code/jetson/requirements-analysis.txt
python -B code/jetson/plot_paper_figure15.py
python -B code/jetson/formal_collection_20261010/plot_power_results.py --run-root data/jetson --temperature
```

The first plotting command rebuilds the paper Figure 15 layout from the verified data (four panels: time, device allocation, module power, and energy per block) in `outputs/jetson/paper/`. Its SVG and PNG are not byte-identical recreations of the typeset PDF. The second command writes separate diagnostic plots to `outputs/jetson/diagnostics/`, including the complete measured-call distribution and the nine raw power traces. Neither command overwrites an existing figure.

## Build and run a new replay

The original platform was an NVIDIA Jetson Orin NX (AArch64, sm_87) in MAXN with CUDA 12.6, driver 540.4.0, and its existing dynamic clock policy. The build requires CMake ≥3.18, a C++17/CUDA17-capable host compiler/toolkit, cuFFT and CUDA Runtime. The formal runner does not use FFTW, Qt, HDF5 or network sockets. Telemetry capture requires the board's `tegrastats` and Linux `/proc`; reading archived records has no such dependency. Nothing here changes the power mode or locks clock frequencies.

Generate the synthetic fixtures with Python, NumPy, SciPy and Matplotlib. This can be done on a host PC. Keep the two generated directories at the documented relative paths in the checkout used for the Jetson build. They are intentionally not tracked; the complete historical fixture bank was about 478 MB, excluding other preflight assets. The recipe uses fixed case definitions and seeds, but regenerated floating-point files and newly built binaries need not be byte-identical across toolchains.

```sh
python -B code/jetson/implementation_v2/python/fixture_export.py --out code/jetson/implementation_v2/python/generated
python -B code/jetson/implementation_v3/python/export_references.py --candidates d4_sinc16

cmake -S code/jetson/formal_collection_20261010/core \
      -B code/jetson/formal_collection_20261010/build \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=87 \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc
cmake --build code/jetson/formal_collection_20261010/build --parallel 2

python3 -B code/jetson/formal_collection_20261010/collect.py freeze
python3 -B code/jetson/formal_collection_20261010/collect.py checks
python3 -B code/jetson/formal_collection_20261010/collect.py latency --run-id replay01
python3 -B code/jetson/formal_collection_20261010/collect.py memory --run-id replay01
python3 -B code/jetson/formal_collection_20261010/collect.py power --run-id replay01
python3 -B code/jetson/formal_collection_20261010/collect.py summarize --run-id replay01
```

The replay controller records the new executable and input hashes before running its 75-case implementation check. It requires that check to pass before acquisition, and does not use the archived permit or executable. Use a fresh checkout/output run ID for a new frozen campaign. Replay build/results directories are git-ignored because they may contain device-local paths and large generated arrays.

Each cached reconstruction begins with a host-ready complex64 block at 400 kHz, 8192 samples/channel (128 Dense channels or 32 CDM/FDM channels). It ends after the complete complex64 output has returned to host memory, with a common 321 × 401 physical grid. The wall timer includes H2D, the common Rx filter, matching, method-specific decimation where applicable, DAS, stitching and D2H completion; file transfer, initialization, CPU output checkpoints, and final destruction are outside that interval. FDM uses fourfold post-match subband decimation and 16-tap local sinc interpolation. Its own-reference check does not claim equal accuracy to the previous full-rate FDM discretization.

Each latency process runs 50 warm-up calls and 200 measured calls; there is no outlier trimming. A power process runs 30 s live-engine idle, 60 s warm-up, then 120 s of submissions followed by the last block's completion. Samples are requested at 100 ms and integrated on actual monotonic timestamps with interpolated boundaries and no idle subtraction. Module-input VDD_IN includes the work of CPU checks and telemetry during the window. The measurement covers numerical-echo reconstruction on the module, rather than acoustic acquisition, a complete sonar frame cycle, or wall-plug power.

## Provenance and scope of edits

The following are historical source-relative locations, not paths required on your machine:

- `jetson_audit_20261010/formal_collection_20261010/core/`, `implementation_v2/core/` and `implementation_v3/core/` under that audit root: source for the measured formal runner. The released CMake, C++/CUDA sources, engine headers, and `telemetry.py` retain their frozen bytes. `contract.json` records the original measurement settings.
- `jetson_audit_20261010/reference_snapshot/numerical_model.py`, and the v2/v3 `python/` directories: deterministic fixture/reference recipes and their original source. The model hash is checked by the exporter.
- `jetson_audit_20261010/formal_collection_20261010/results/main_20261010/`: raw logs and summaries, now organized under `data/jetson/`. The numeric logs are unmodified. Only the nine invocation records were reduced to public phase events and status, omitting machine-local command paths, stdout and process IDs.
- `RINENG_review/latex/standalone_v08_20261014/main/main.tex`, the `tab:runtime` and `fig10` labels: manuscript table/figure provenance. The historic Figure 15 asset is byte-identical to the v01 `figure10.pdf`, whose plotting recipe is `RINENG_review/work/v01_20261010/plot_device_pipeline.py`. The released four-panel recipe reads the measured memory ledger instead of hardcoding its figures.

The release controller has a portability adaptation to generate new manifests/receipts without private historical audit inputs. The diagnostic plot script has portable output/provenance paths. Those two scripts are not presented as byte-identical to the original acquisition controller/plotter. Historical manifests are kept only as provenance; they intentionally name historical assets not redistributed here. This release contains no original GUI project, toolchain binary, SSH helper, credential, host-key database, or generated fixture archive.

Package QA: offline table verification passed; both fixture exporters and the replay controller were checked for import/CLI availability. This checkout was not rebuilt or executed on Jetson during release preparation.
