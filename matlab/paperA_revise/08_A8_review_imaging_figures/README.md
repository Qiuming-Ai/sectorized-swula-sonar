# A8 historical controlled imaging export

A8 is an adapter around six historical Exp3/Exp4 scripts, not a standalone
beamformer. It reads and applies the original controlled text transformations
to the source in memory, then captures images at the selected SNRs. Its
single-trial, per-method measured AWGN and archived geometry are not the final
paper's common per-element noise and 30-trial image results.

Run from the repository root:

```sh
matlab -batch "addpath('matlab/paperA_revise/08_A8_review_imaging_figures'); export_review_imaging_data()"
python matlab/paperA_revise/08_A8_review_imaging_figures/plot_review_imaging.py
python matlab/paperA_revise/08_A8_review_imaging_figures/compute_review_image_metrics.py
```

The MATLAB entry uses the sibling `matlab/legacy_sonarsim/` archive by
default. An explicit trusted compatible directory can be passed instead:
`export_review_imaging_data('/path/to/SonarSim')`. Preflight checks the six
required entry scripts before creating any output. The adapter evaluates
MATLAB source, so only pass an author-provided or otherwise trusted source
tree. Historical source text must retain the patterns used for data capture
and the stated controlled settings.

The six scripts are the Dense, CDM and FDM
`PaperA_exp3_*_noise_doppler_baseband.m` entry points and the random,
Nested and Coprime `PaperA_exp4_*_noise_doppler_baseband.m` entry points.
They require the matching archived helper functions and toolbox functions
such as Signal Processing Toolbox's `hilbert`/`resample` and Communications
Toolbox's `awgn`. The CDM branch calls the legacy
`sim_rx_from_scatterers_perTX_cuda` path and requires its supported runtime
or compiled MEX as specified in the legacy archive. See that archive's
README for asset/build boundaries. A8 cannot run a complete matched
comparison if one of those dependencies is unavailable.

The entry writes data to this folder's `data/` and transformed diagnostic
scripts to `debug_patched/`. Python reads these generated data; no historical
MAT, PNG, PDF, or generated debug script is bundled.
