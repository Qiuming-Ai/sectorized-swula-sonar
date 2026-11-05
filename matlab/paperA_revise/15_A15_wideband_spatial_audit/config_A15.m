function cfg=config_A15(mode)
%CONFIG_A15 Centre-frequency versus waveform-weighted spatial response.
if nargin<1,mode="quick";end
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'common'));
cfg=exp_default_config(mode);
if lower(string(mode))=="quick"
    cfg.angle_grid_deg=-50:0.1:50;cfg.frequency_samples=101;
else
    cfg.angle_grid_deg=-50:0.025:50;cfg.frequency_samples=401;
end
cfg.target_angles_deg=[-41.25,-24.25,-8,8,24.25,41.25];
cfg.methods=["cdm","fdm_current"];
cfg.method_labels=["CDM--ZC","packed FDM--LFM"];
cfg.frequency_grid_hz=linspace(cfg.fc-cfg.total_bw/2,cfg.fc+cfg.total_bw/2,cfg.frequency_samples);
cfg.reference_model="waveform-power-weighted coherent integration of exact per-frequency transmit and receive propagation phases with TTD steering";
cfg.approximation_model="single spatial manifold at each waveform's acoustic centre frequency";
end
