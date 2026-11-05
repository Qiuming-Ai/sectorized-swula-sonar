function cfg = exp_default_config(mode)
%EXP_DEFAULT_CONFIG Shared physical and numerical settings for Paper A revisions.
%   cfg = EXP_DEFAULT_CONFIG("quick") is intended for code checks.
%   cfg = EXP_DEFAULT_CONFIG("paper") uses finer grids and 30 Monte-Carlo runs.

if nargin < 1
    mode = "quick";
end
mode = lower(string(mode));
if ~ismember(mode, ["quick", "paper"])
    error('mode must be "quick" or "paper".');
end

cfg.mode = mode;
cfg.c = 1500;
cfg.fc = 450e3;
cfg.total_bw = 100e3;
cfg.frequency_grid = (400:25:500) * 1e3;

cfg.n_tx = 12;
cfg.n_rx_sparse = 32;
cfg.tx_spacing_lambda = 0.5;
cfg.rx_spacing_lambda = 2.0;
cfg.n_rx_dense = 1 + round((cfg.n_rx_sparse - 1) * ...
    cfg.rx_spacing_lambda / 0.5); % same aperture, lambda/2 sampling
cfg.dense_rx_spacing_lambda = 0.5;

% Revised 100-degree FOV. The six steering angles were optimized jointly
% over 400--500 kHz for worst-case coverage with the existing 12-element Tx.
cfg.fov_deg = [-50, 50];
cfg.focus_angles_deg = [-41.25, -24.25, -8.00, 8.00, 24.25, 41.25];
sector_breaks = linspace(cfg.fov_deg(1),cfg.fov_deg(2),7);
cfg.sector_edges_deg = [sector_breaks(1:end-1).',sector_breaks(2:end).'];

% Fixed physical dimensions, expressed relative to lambda at 450 kHz.
% Replace these two assumptions with measured element widths if available.
cfg.tx_element_width_lambda = 0.45;
cfg.rx_element_width_lambda = 0.45;
cfg.include_element_pattern = true;
cfg.steering_mode = "ttd"; % "ttd" or fixed phase shifters at cfg.fc

% Original-paper waveforms. The manuscript uses a 2-MHz acquisition rate;
% this reduced-order simulator retains a 400-kHz complex-baseband rate,
% which is sufficient for the 100-kHz processed band.
cfg.acquisition_fs_hz = 2e6;
cfg.waveform_fs = 400e3;
cfg.pulse_duration_s = 5.0e-3;
% L=499 is the nearest prime-length ZC code to B*T=500 chips. It gives a
% 4.99-ms coded pulse at the 100-kHz chip rate while preserving CAZAC roots.
cfg.zc_length = 499;
cfg.zc_roots = [5, 7, 1, 3, 11, 9];
cfg.zc_chip_rate_hz = cfg.total_bw;
cfg.waveform_definition = "original_paper_fdm_lfm_and_cdm_zc";
cfg.max_range_lag_samples = 64;
cfg.current_fdm_centers_hz = [426, 442, 490, 474, 458, 412] * 1e3;
cfg.current_fdm_bandwidth_hz = 20e3;
cfg.fdm_guard_fraction = 0.10;
cfg.fdm_optimized_permutation = 1:6;

if mode == "quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):1.0:cfg.fov_deg(2);
    cfg.n_range_samples = 161;
    cfg.monte_carlo_runs = 5;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.25:cfg.fov_deg(2);
    cfg.n_range_samples = 321;
    cfg.monte_carlo_runs = 30;
end
cfg.range_sample_spacing_m = cfg.c / (2 * cfg.waveform_fs);
cfg.range_axis_m = (0:cfg.n_range_samples-1) * cfg.range_sample_spacing_m;
cfg.metric_guard_radius = [2, 1]; % [range samples, angle samples]
cfg.default_seed = 20260810;
cfg.default_snr_db = 20;
end
