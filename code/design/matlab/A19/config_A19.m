function cfg=config_A19(mode)
%CONFIG_A19 Broadband sector-count coverage and alias constraints.
if nargin<1,mode="quick";end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg=exp_default_config(mode);
cfg.candidate_sector_counts=4:8;
cfg.coverage_threshold_db=-6.6;
cfg.grating_threshold_db=-24.0;
if lower(string(mode))=="quick"
    cfg.design_angle_step_deg=0.25;
    cfg.search_audit_angle_step_deg=0.50;
    cfg.search_frequency_grid=(400:10:500)*1e3;
    cfg.final_audit_angle_step_deg=0.10;
    cfg.final_audit_frequency_grid=(400:5:500)*1e3;
else
    cfg.design_angle_step_deg=0.05;
    cfg.search_audit_angle_step_deg=0.25;
    cfg.search_frequency_grid=(400:5:500)*1e3;
    cfg.final_audit_angle_step_deg=0.01;
    cfg.final_audit_frequency_grid=(400:1:500)*1e3;
end
cfg.response_angle_grid_deg=-90:cfg.search_audit_angle_step_deg:90;
cfg.report_angle_grid_deg=cfg.fov_deg(1):cfg.search_audit_angle_step_deg:cfg.fov_deg(2);
cfg.receive_alias_orders=-3:3;
cfg.receive_alias_orders(cfg.receive_alias_orders==0)=[];
end
