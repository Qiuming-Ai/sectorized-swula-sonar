function cfg = config_A1(mode)
%CONFIG_A1 Configuration for the exact transmit beampattern experiment.
if nargin < 1, mode = "quick"; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
if cfg.mode == "paper"
    cfg.angle_grid_deg = -90:0.05:90;
    cfg.frequency_grid = linspace(400e3,500e3,21);
else
    cfg.angle_grid_deg = -90:0.10:90;
end
end
