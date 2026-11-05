function cfg=config_A14(mode)
%CONFIG_A14 Pure receive versus effective two-way response audit.
if nargin<1,mode="quick";end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'common'));
cfg=exp_default_config(mode);
if lower(string(mode))=="quick",cfg.angle_grid_deg=-60:0.1:60;else,cfg.angle_grid_deg=-60:0.02:60;end
cfg.target_angles_deg=[0,8,41.25];
cfg.methods=["dense","cdm","fdm_current"];
cfg.method_labels=["Dense ULA","CDM SW-ULA","packed FDM SW-ULA"];
cfg.mainlobe_definition="contiguous lobe bounded by the first local minima around the peak";
end
