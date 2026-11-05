function cfg = config_A13(mode)
%CONFIG_A13 Unified raw-channel Monte Carlo metric protocol.
if nargin<1, mode="quick"; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg=exp_default_config(mode);
cfg.methods=["dense","cdm","fdm_current"];
cfg.method_labels=["Dense ULA","CDM SW-ULA","packed FDM SW-ULA"];
cfg.scenario_names=["point18","dense1000x50"];
cfg.snr_values_db=[-15,-10,-5,0];
cfg.n_rx_dense=128;
cfg.angle_grid_deg=cfg.fov_deg(1):1:cfg.fov_deg(2);
cfg.n_range_samples=321;
cfg.range_axis_m=(0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
cfg.range_origin_m=4.92;
cfg.physical_range_axis_m=cfg.range_origin_m+cfg.range_axis_m;
cfg.metric_core_radius=[2,1];
cfg.metric_guard_radius=[5,2];
cfg.scene_seed=2026081201;
cfg.noise_seed_base=2026081300;
if lower(string(mode))=="quick", cfg.monte_carlo_runs=3; else, cfg.monte_carlo_runs=30; end
cfg.snr_definition="mean noise-free raw complex-baseband channel power / mean added channel-noise power";
cfg.pairing_rule="shared standardized noise samples on common sample/channel indices, with per-method scaling to the declared aggregate input SNR";
cfg.cnr_scope="reported for the 18-point scene; descriptive only for dense support because 50 range layers leave no independent empty-water background model";
cfg.nmse_scope="complex noise-induced NMSE relative to the same method's noise-free reconstruction after one complex scalar fit";
cfg.dense_reference_correlation_scope="linear-magnitude Pearson correlation to the paired dense noise-free reconstruction of the same scene";
end
