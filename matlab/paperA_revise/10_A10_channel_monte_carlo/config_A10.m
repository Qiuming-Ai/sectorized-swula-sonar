function cfg = config_A10(mode)
%CONFIG_A10 Raw complex-baseband paired-noise Monte-Carlo protocol.
if nargin<1, mode="quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.methods = ["dense","cdm","fdm_current"];
cfg.method_labels = ["Dense ULA","CDM SW-ULA","FDM SW-ULA"];
cfg.scenario_names = ["point18","dense1000x50"];
cfg.snr_values_db = [-15,-10,-5,0];
cfg.n_rx_dense = 128;
cfg.angle_grid_deg = cfg.fov_deg(1):1:cfg.fov_deg(2);
cfg.n_range_samples = 321;
cfg.range_axis_m = (0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
cfg.range_origin_m = 4.92;
cfg.physical_range_axis_m = cfg.range_origin_m+cfg.range_axis_m;
cfg.metric_guard_radius = [5,2];
cfg.scene_seed = 2026081201;
cfg.noise_seed_base = 2026081300;
if lower(string(mode))=="quick"
    cfg.monte_carlo_runs = 3;
else
    cfg.monte_carlo_runs = 30;
end
cfg.snr_definition = "mean noise-free raw complex-baseband channel power / mean added channel-noise power";
cfg.pairing_rule = "same standardized complex-noise realization for common samples/channels; per-method scaling to equal input SNR";
end
