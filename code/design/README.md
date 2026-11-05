# Array and waveform design figures

This directory contains the design examples and historical single-frequency
diagnostics for Figures 1–8. It is separate from the implemented common-FIR
waveforms, full two-way response, and imaging experiments. A plot from the
distributed numerical data need not reproduce the manuscript PDF layout or
file hash.

## Quick redraw

From the repository root, with Python, NumPy, SciPy, and Matplotlib installed:

```sh
python code/design/plot_design.py --verify-only
python code/design/plot_design.py --output outputs/design
python code/design/plot_design.py --figures 3 8 --output outputs/design
```

These commands use the recorded curve data or the stated analytical
expressions and do not access the original workspace, a MATLAB installation,
or a manuscript PDF.

| Figure | What can be reproduced |
|---|---|
| 1: receive-array aliasing | Historical raster composite. `PaperA_exp0.m` is a related one-spacing 5.12-m angular-cut script with the matching target directions, but the three underlying historical traces/settings were not retained together in a verified source set. Not included as a runnable exact figure. |
| 2: transmit masking | Numerical slices from `PSF.mat` using the zero-based equivalents of rows 851, 791, and 911 of the MATLAB recipe `PSF.m`. `data/design/figure2_legacy_psf_slices.csv` retains the linear normalized values before the 20-log display conversion. The old annotations were assembled separately; the raw MAT file and the exact composite are not distributed here. |
| 3: six transmit beams | Archived final-figure curve samples in `data/spatial/Figure03_transmit_beampattern_data.csv`, and A1 configuration/metrics in `data/design`. The A1 MATLAB scripts independently generate the 450-kHz beampatterns and their original frequency sweep. |
| 4(a): FDM design | Zero-motion **prototype** in the legacy `window1.m`, plus a Python redraw of its correlation. It uses a 0.5-ms pulse, 20-kHz chirp, and the archived common cosine envelope multiplied by the indicated additional taper. `Rect` denotes the additional rectangular multiplier, not removal of the common envelope. The legacy variable named `vel_list` indexes window choices here: the echo calls themselves use zero velocity. This example is not the 5-ms, FIR-filtered implementation in Fig. 5. |
| 4(b): ideal CDM design | Historical `CDM_ZC_crosscorr_ratio.m` and a Python implementation of the same ideal periodic amplitude bound, `20 log10(M−1) − 10 log10(L)`. It illustrates scaling under the appropriate distinct-root correlation conditions; it is not an aperiodic finite-pulse evaluation or a guarantee for arbitrary composite lengths. The original script evaluates longer durations than the visible 0–10 ms redraw. |
| 5: finite-pulse ambiguity | See the separate waveform/ambiguity reproducibility entry. Do not use Fig. 4 prototype curves as substitutes. |
| 6: simultaneous multibeam architecture | Historical raster schematic; no verified editable source was located. No scientific curve data are implied and no guessed rebuild is included. |
| 7: DAS pipeline | `data/design/figure7_pipeline_geometry.json` preserves the eight verified stage labels and 22 vector paths (eight nodes plus seven arrow stems and heads). The portable redraw uses these paths and a local serif font. This diagram is not the executable reconstruction chain. |
| 8: angular spatial product | `data/spatial/A14_angular_profiles.csv` contains the original 0.02° angular samples, and `data/design/A14_angular_metrics.csv` is the source of the single-frequency spatial-product metrics in manuscript Table 2. The plot selects the five 8° traces and displays only ±50°. A14 is a historical matched-aperture diagnostic (125 Dense elements in its recorded configuration); the later common-FIR scene simulations use the 128-element Dense reference in Table 1. The archived A14 narrowband/gated profiles are not the full simultaneous-waveform-bank response in Fig. 9. |

The Fig. 2 and Fig. 7 extraction receipts contain hashes of the archived source
files, not copies of the original manuscript graphics. Re-extraction is
optional for holders of those archives and additionally requires h5py and
PyMuPDF:

```sh
python code/design/export_design_inputs.py --legacy-mat path/to/PSF.mat --pipeline-pdf path/to/original_pipeline.pdf --output outputs/design-inputs
```

## Historical MATLAB design audits

The archived A1 and A19 scripts and only their required helpers are under
`matlab/`; they are included unchanged. With MATLAB R2025b, run from the
repository root:

```sh
matlab -batch "addpath('code/design/matlab/A1'); run_A1_tx_beampattern('paper')"
matlab -batch "addpath('code/design/matlab/A19'); run_A19_sector_count_optimization('paper')"
```

Both accept `'quick'` for a smaller exploratory run. Generated results are
placed in a new `results/` subfolder beside the corresponding entry point.
Frozen results in `data/design` are not overwritten.

A19 is the historical **per-sector** steering search and independent audit of
the implemented six-sector settings; it is not the later joint constrained
selection. In the archived A19 `paper` code, the candidate-scoring loop uses
`cfg.frequency_grid` (400:25:500 kHz), even though the exported
`SearchFrequencySpacingKHz` field describes the configured
`search_frequency_grid` (400:5:500 kHz). This mismatch is retained rather than
silently changing the archived experiment. The final audit separately uses its
declared 0.01° / 1-kHz grid with inserted sector boundaries. The implemented
six-sector row is −6.51734 dB worst assigned coverage and −24.32507 dB worst
visible-alias ratio. Use the separate joint coverage search for the final
manuscript's search protocol and candidate claims.

Table 1's actual pulse and scene parameters are defined by the release's
common numerical model, not by the historical A1/A14/A19 configurations.
Table 3's effective-template width and sidelobe metrics come from
`data/waveforms/window/window_template_metrics.csv`.
