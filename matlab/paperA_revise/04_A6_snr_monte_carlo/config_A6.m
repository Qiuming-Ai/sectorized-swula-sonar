function cfg = config_A6(mode)
%CONFIG_A6 SNR and Monte-Carlo settings.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.snr_values_db = [-20,-15,-10,-5,0,5,10];
cfg.methods = ["dense","cdm","fdm_current","fdm_strict"];
cfg.scenario_names = ["sparse","medium","dense"];
cfg.scenario_scatter_counts = [10,200,1000];
cfg.scenario_range_layers = [2,10,50];
cfg.scenario_dynamic_ranges_db = [10,20,40];
if cfg.mode == "quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):2:cfg.fov_deg(2);
    cfg.n_range_samples = 121;
    cfg.monte_carlo_runs = 5;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.5:cfg.fov_deg(2);
    cfg.n_range_samples = 241;
    cfg.monte_carlo_runs = 30;
end
cfg.range_axis_m = (0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
cfg.max_range_lag_samples = cfg.n_range_samples-1;
end
