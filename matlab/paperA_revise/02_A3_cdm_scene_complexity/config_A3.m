function cfg = config_A3(mode)
%CONFIG_A3 Scene-complexity sweep settings.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.scatter_counts = [1,10,50,200,500,1000];
cfg.range_layer_counts = [1,5,10,20,50];
cfg.dynamic_ranges_db = [0,10,20,40];
cfg.methods = ["dense","cdm","fdm_current","fdm_strict"];
cfg.reference_scatter_count = 200;
cfg.reference_layer_count = 10;
cfg.reference_dynamic_range_db = 20;
cfg.full_factorial = false; % true is much slower; OFAT covers every requested level
if cfg.mode == "quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):2:cfg.fov_deg(2);
    cfg.n_range_samples = 121;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.5:cfg.fov_deg(2);
    cfg.n_range_samples = 241;
end
cfg.range_axis_m = (0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
% Use the complete scene-delay window; short-lag checks can hide the
% accumulation of range sidelobes in a 50-layer diffuse scene.
cfg.max_range_lag_samples = cfg.n_range_samples-1;
end
