function cfg = config_A4(mode)
%CONFIG_A4 FDM band-layout comparison settings.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.methods = ["fdm_current","fdm_strict","fdm_optimized"];
if cfg.mode == "quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):1:cfg.fov_deg(2);
    cfg.n_range_samples = 161;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.25:cfg.fov_deg(2);
    cfg.n_range_samples = 321;
end
cfg.range_axis_m = (0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
cfg.max_range_lag_samples = cfg.n_range_samples-1;
cfg.optimize_assignment = true;
end
