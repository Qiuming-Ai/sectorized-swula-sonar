function cfg = config_A7(mode)
%CONFIG_A7 Fair array-baseline settings.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
if cfg.mode == "quick"
    cfg.angle_grid_deg = -60:0.5:60;
    cfg.frequency_grid = (400:25:500)*1e3;
else
    cfg.angle_grid_deg = -60:0.1:60;
    cfg.frequency_grid = linspace(400e3,500e3,21);
end
cfg.point_target_angles_deg = [-40,-15,0,20,45];
cfg.adaptive_source_angles_deg = [-30,5,38];
cfg.adaptive_source_powers = [1,0.7,0.5];
cfg.adaptive_snr_db = 10;
cfg.diffuse_scatter_count = 100;
end
