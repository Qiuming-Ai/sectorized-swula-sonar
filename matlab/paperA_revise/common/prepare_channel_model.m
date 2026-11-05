function model = prepare_channel_model(cfg, method)
%PREPARE_CHANNEL_MODEL Element-level complex-baseband forward/receive model.
%   The model retains raw receive channels until after additive noise is
%   injected.  It is a deterministic 2-D far-field channel simulation, not
%   a propagation, reverberation, or measured-transducer model.

method = lower(string(method));
angles = cfg.angle_grid_deg(:).';
lambda_c = cfg.c/cfg.fc;

switch method
    case "dense"
        n_rx = cfg.n_rx_dense;
        spacing_lambda = cfg.dense_rx_spacing_lambda;
        [bank, meta] = make_waveform_bank(cfg,"dense");
        waveforms = bank(:,1);
        center_frequency_hz = cfg.fc;
        tx_response = {ones(1,numel(angles))};
        output_sector = {true(1,numel(angles))};
    case {"cdm","fdm_current","fdm_strict","fdm_optimized"}
        n_rx = cfg.n_rx_sparse;
        spacing_lambda = cfg.rx_spacing_lambda;
        [waveforms, meta] = make_waveform_bank(cfg,method);
        % Equalize total waveform energy across the simultaneous six-beam
        % transmission and the single-waveform dense reference.
        waveforms = waveforms/sqrt(size(waveforms,2));
        center_frequency_hz = meta.center_frequency_hz;
        masks = sector_masks(cfg,angles);
        tx_response = cell(1,size(waveforms,2));
        output_sector = cell(1,size(waveforms,2));
        for k = 1:size(waveforms,2)
            [tx,~] = array_responses(cfg,center_frequency_hz(k),angles);
            tx_response{k} = tx(k,:);
            output_sector{k} = masks(k,:);
        end
    otherwise
        error('Unknown method: %s',method);
end

rx_positions_m = ((0:n_rx-1)-(n_rx-1)/2).' * spacing_lambda*lambda_c;
rx_taper = local_window("hamming",n_rx);
rx_taper = rx_taper/sum(rx_taper);

n_wave = size(waveforms,2);
rx_manifold = cell(1,n_wave);
rx_weights = cell(1,n_wave);
for k = 1:n_wave
    frequency = center_frequency_hz(min(k,numel(center_frequency_hz)));
    wavenumber = 2*pi*frequency/cfg.c;
    direction_cosine = sind(angles);
    manifold = exp(1j*wavenumber*rx_positions_m*direction_cosine);
    if isfield(cfg,'include_element_pattern') && cfg.include_element_pattern
        width_m = cfg.rx_element_width_lambda*lambda_c;
        x = (width_m*frequency/cfg.c).*direction_cosine;
        element = ones(size(x));
        nz = abs(x)>eps;
        element(nz) = sin(pi*x(nz))./(pi*x(nz));
        manifold = manifold.*element;
    end
    rx_manifold{k} = manifold;
    rx_weights{k} = rx_taper.*exp(-1j*wavenumber*rx_positions_m*direction_cosine);
end

model.method = method;
model.waveforms = waveforms;
model.waveform_meta = meta;
model.tx_response = tx_response;
model.output_sector = output_sector;
model.rx_manifold = rx_manifold;
model.rx_weights = rx_weights;
model.rx_positions_m = rx_positions_m;
model.n_rx = n_rx;
model.n_waveforms = n_wave;
model.n_range = cfg.n_range_samples;
model.n_angle = numel(angles);
model.noise_location = "complex receive channels before matched filtering and beamforming";
model.energy_normalization = "unit total discrete waveform energy";
end
