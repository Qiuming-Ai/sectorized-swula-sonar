function model = prepare_reduced_model(cfg, method)
%PREPARE_REDUCED_MODEL Precompute invariant waveform and spatial operators.
method = lower(string(method));
angles = cfg.angle_grid_deg;
m = numel(cfg.focus_angles_deg);
switch method
    case "dense", waveform_family = "dense";
    case "cdm", waveform_family = "cdm";
    case {"fdm_strict","fdm_current","fdm_optimized"}
        waveform_family = method;
    otherwise
        error('Unknown method: %s',method);
end

[waveforms,wave_meta] = make_waveform_bank(cfg,waveform_family);
[kernels,cross_peak,lags] = xcorr_kernel_bank(waveforms,cfg.max_range_lag_samples);
model.method = method;
model.waveform = wave_meta;
model.kernels = kernels;
model.cross_peak = cross_peak;
model.lags = lags;
model.masks = sector_masks(cfg,angles);
model.source_sets = cell(1,m);
for k = 1:m
    model.source_sets{k} = find(wave_meta.slot_index==wave_meta.slot_index(k));
end

if method == "dense"
    [~,model.dense_rx] = array_responses(cfg,cfg.fc,angles, ...
        Nrx=cfg.n_rx_dense,RxSpacingLambda=cfg.dense_rx_spacing_lambda);
    return;
end

tx_by_source = complex(zeros(m,numel(angles)));
rx_by_output = cell(1,m);
for k = 1:m
    frequency = wave_meta.center_frequency_hz(k);
    [tx,rx] = array_responses(cfg,frequency,angles);
    tx_by_source(k,:) = tx(k,:);
    rx_by_output{k} = rx;
end
model.spatial = cell(m,m);
for k = 1:m
    scan_idx = model.masks(k,:);
    for i = 1:m
        model.spatial{k,i} = rx_by_output{k}(scan_idx,:) .* tx_by_source(i,:);
    end
end
end
