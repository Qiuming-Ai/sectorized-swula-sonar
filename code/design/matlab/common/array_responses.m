function [tx_response, rx_response, meta] = array_responses(cfg, frequency_hz, angles_deg, options)
%ARRAY_RESPONSES Complex transmit and receive responses for fixed ULAs.
% tx_response is Nbeam-by-Nangle; rx_response is Nscan-by-Ntarget.

arguments
    cfg struct
    frequency_hz (1,1) double {mustBePositive}
    angles_deg (1,:) double
    options.Nrx (1,1) double {mustBeInteger,mustBePositive} = cfg.n_rx_sparse
    options.RxSpacingLambda (1,1) double {mustBePositive} = cfg.rx_spacing_lambda
    options.TxWindow = "hamming"
    options.RxWindow = "hamming"
end

c = cfg.c;
lambda_c = c / cfg.fc;
k = 2*pi*frequency_hz/c;
if lower(string(cfg.steering_mode)) == "ttd"
    k_steer = k;
else
    k_steer = 2*pi*cfg.fc/c;
end

ptx = ((0:cfg.n_tx-1) - (cfg.n_tx-1)/2).' * ...
    cfg.tx_spacing_lambda * lambda_c;
prx = ((0:options.Nrx-1) - (options.Nrx-1)/2).' * ...
    options.RxSpacingLambda * lambda_c;
wtx = local_window(options.TxWindow, cfg.n_tx);
wrx = local_window(options.RxWindow, options.Nrx);
wtx = wtx / sum(wtx);
wrx = wrx / sum(wrx);

u = sind(angles_deg);
tx_response = complex(zeros(numel(cfg.focus_angles_deg), numel(angles_deg)));
for ib = 1:numel(cfg.focus_angles_deg)
    steer = exp(-1j*k_steer*ptx*sind(cfg.focus_angles_deg(ib)));
    tx_response(ib,:) = (wtx .* steer).' * exp(1j*k*ptx*u);
end

rx_response = complex(zeros(numel(angles_deg), numel(angles_deg)));
target_steering = exp(1j*k*prx*u);
for iscan = 1:numel(angles_deg)
    scan_weight = wrx .* exp(-1j*k_steer*prx*u(iscan));
    rx_response(iscan,:) = scan_weight.' * target_steering;
end

if isfield(cfg, 'include_element_pattern') && cfg.include_element_pattern
    tx_width_m = cfg.tx_element_width_lambda * lambda_c;
    rx_width_m = cfg.rx_element_width_lambda * lambda_c;
    tx_element = normalized_sinc((tx_width_m*frequency_hz/c) .* u);
    rx_element = normalized_sinc((rx_width_m*frequency_hz/c) .* u);
    tx_response = tx_response .* tx_element;
    rx_response = rx_response .* rx_element;
else
    tx_element = ones(size(u));
    rx_element = ones(size(u));
end

meta.frequency_hz = frequency_hz;
meta.tx_positions_m = ptx;
meta.rx_positions_m = prx;
meta.tx_element_pattern = tx_element;
meta.rx_element_pattern = rx_element;
end

function y = normalized_sinc(x)
y = ones(size(x));
nz = abs(x) > eps;
y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
end
