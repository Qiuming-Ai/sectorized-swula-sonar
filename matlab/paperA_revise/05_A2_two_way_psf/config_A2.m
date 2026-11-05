function cfg = config_A2(mode)
%CONFIG_A2 Two-way broadband point-spread settings.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.psf_methods = ["dense_lambda2","sw_ula_wideband", ...
    "packed_fdm","sector_fdm","cdm_zc","ideal_cdm"];
cfg.psf_target_angles_deg = [8.0,41.25];
cfg.psf_guard = [5,2]; % [range-lag samples, angle samples]
cfg.max_range_lag_samples = 240;
if cfg.mode == "quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):0.5:cfg.fov_deg(2);
    cfg.frequency_grid = (400:25:500)*1e3;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.1:cfg.fov_deg(2);
    cfg.frequency_grid = linspace(400e3,500e3,21);
end
end
