function cfg = config_A5(mode)
%CONFIG_A5 Doppler-sensitivity settings.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.velocities_mps = [-2,-1,-0.5,0,0.5,1,2];
cfg.profile_velocities_mps = [0,0.5,1];
cfg.profile_beam = 3;
cfg.methods = ["cdm","fdm_current","fdm_strict"];
cfg.doppler_model = "wideband_time_scale_plus_residual_carrier";
if cfg.mode == "quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):2:cfg.fov_deg(2);
    cfg.n_range_samples = 121;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.5:cfg.fov_deg(2);
    cfg.n_range_samples = 241;
end
cfg.range_axis_m = (0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
cfg.max_range_lag_samples = cfg.n_range_samples-1;
end
