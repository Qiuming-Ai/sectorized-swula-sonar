# MATLAB figure prototypes

These are study-related scripts from the historical `paperA` figure workspace, not the final waveform or imaging validation.

- `PSF.m` redraws the Figure 2 explanatory curves from the shipped normalized slices. The optional original `PSF.mat` workspace can still be used if present beside the script.
- `window1.m` and `CDM_ZC_crosscorr_ratio.m` generate the historical Figure 4 waveform-design examples. `doppler1.m` and `lfm_pulse_compression.m` are early waveform illustrations, not substitutes for the final Figure 5 data.
- `draw_farfield_active_point_scat.m` and `plot_scan_fov_vs_spacing.m` are explanatory geometry/array scripts. `SignalFilter.m` and `fd_fir_delay_selftest.m` retain the early filter implementation and its fractional-delay example.

Run from a working directory where generated figures may be written. Add only this folder to the MATLAB path for these examples; there are similarly named historical scripts in other modules. Proprietary toolboxes and unrelated File Exchange graphics demos are not redistributed.
