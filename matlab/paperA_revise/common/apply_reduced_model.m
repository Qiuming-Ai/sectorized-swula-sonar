function [image_complex,diagnostics] = apply_reduced_model(cfg,scene,model,snr_db,seed)
%APPLY_REDUCED_MODEL Apply a model prepared by PREPARE_REDUCED_MODEL.
if nargin < 4 || isempty(snr_db), snr_db = Inf; end
if nargin < 5, seed = cfg.default_seed; end
nr = cfg.n_range_samples;
na = numel(cfg.angle_grid_deg);
if ~isequal(size(scene),[nr,na])
    error('scene size is inconsistent with cfg.');
end

if model.method == "dense"
    projected = scene * model.dense_rx.';
    image_complex = conv2(projected,model.kernels(:,1,1),'same');
else
    m = numel(cfg.focus_angles_deg);
    image_complex = complex(zeros(nr,na));
    for k = 1:m
        scan_idx = model.masks(k,:);
        accum = complex(zeros(nr,nnz(scan_idx)));
        source_set = model.source_sets{k};
        for i = source_set
            projected = scene * model.spatial{k,i}.';
            accum = accum + conv2(projected,model.kernels(:,i,k),'same');
        end
        image_complex(:,scan_idx) = accum;
    end
end

signal_power = mean(abs(image_complex(:)).^2);
if isfinite(snr_db) && signal_power > 0
    rng(seed,'twister');
    noise_power = signal_power/10^(snr_db/10);
    image_complex = image_complex + sqrt(noise_power/2).* ...
        (randn(size(image_complex))+1j*randn(size(image_complex)));
else
    noise_power = 0;
end
diagnostics.method = model.method;
diagnostics.waveform = model.waveform;
diagnostics.cross_peak = model.cross_peak;
diagnostics.lags = model.lags;
diagnostics.output_signal_power = signal_power;
diagnostics.output_noise_power = noise_power;
diagnostics.snr_definition = "mean reconstructed-output signal power / noise power";
end
