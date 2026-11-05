# A16 — Shipwreck scattering-phase ensemble

This is the historical MATLAB phase-ensemble control, **not** the final v08 Python Figure 14 or its confidence intervals. It uses independent uniform unresolved-scatterer phases on the same one-bounce arrival catalogue, paired standardized noise, the historical equal-input-SNR channel chain, and correlation to the paired noisy Dense image. The legacy interval calculation is `1.96*SD/sqrt(N)` and the guards are in grid bins.

From the package root in MATLAB:

```matlab
addpath(fullfile(pwd,'matlab','paperA_revise','16_A16_shipwreck_phase_ensemble'));
summary = run_A16_shipwreck_phase_ensemble("quick");  % 3 phases, full archive grid
```

`"paper"` selects 30 phases of the historical MATLAB chain. By default A16 rebuilds the common legacy A11 scene from the archived arrival CSV and STL geometry catalogue; it does not retrace the STL or require a prebuilt workspace/MEX. The optional arguments are `geometry_dir`, `out_dir` and `saved_a11_workspace`. Use the last one only for a compatible **paper-mode** A11 workspace. Default outputs are `outputs/matlab/A16_quick/` or `outputs/matlab/A16_paper/`; the diagnostic plotter accepts `--results PATH_TO_RUN`. The port records each method's own SNR metadata rather than copying the last method's diagnostic into every row; image/metric computations are unchanged.

For the final paper figures and 30-repetition, Student-t reference statistics use the separate released Python scripts and frozen trial data. The historical MATLAB controls are numerical simulations, not acoustic acquisition.
