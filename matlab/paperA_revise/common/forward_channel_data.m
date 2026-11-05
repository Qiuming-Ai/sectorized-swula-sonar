function [raw_clean, diagnostics] = forward_channel_data(cfg, scene, model)
%FORWARD_CHANNEL_DATA Generate noise-free element-level baseband channels.

if ~isequal(size(scene),[cfg.n_range_samples,numel(cfg.angle_grid_deg)])
    error('scene size is inconsistent with cfg.');
end

n_wave_sample = size(model.waveforms,1);
n_raw = cfg.n_range_samples+n_wave_sample-1;
raw_clean = complex(zeros(n_raw,model.n_rx));

for k = 1:model.n_waveforms
    angle_weighted_scene = scene.*model.tx_response{k};
    element_impulses = angle_weighted_scene*model.rx_manifold{k}.';
    raw_clean = raw_clean + fft_filter_channels(element_impulses, ...
        model.waveforms(:,k),n_raw);
end

diagnostics.method = model.method;
diagnostics.raw_signal_power = mean(abs(raw_clean(:)).^2);
diagnostics.raw_peak = max(abs(raw_clean(:)));
diagnostics.n_raw_samples = size(raw_clean,1);
diagnostics.n_receive_channels = size(raw_clean,2);
diagnostics.noise_location = model.noise_location;
diagnostics.energy_normalization = model.energy_normalization;
end
