function [kernels,peak_matrix,peak_lag_samples,lags] = subcarrier_kernel_bank( ...
    cfg,waveforms,velocity_mps,max_lag,compensation_velocity_mps)
%SUBCARRIER_KERNEL_BANK Circular spectral-demultiplexing ambiguity kernels.
% At zero Doppler, disjoint DFT supports have zero cross terms for every
% circular delay. Doppler time scaling causes energy leakage across supports.
% An optional matched velocity hypothesis removes this scale and carrier
% offset within the same received frame; it is not slow-time/TDM processing.
if nargin < 5, compensation_velocity_mps = 0; end
[n,m] = size(waveforms);
max_lag = min(max_lag,floor(n/2)-1);
lags = (-max_lag:max_lag).';
t = (0:n-1).'/cfg.waveform_fs;
eta = (cfg.c+velocity_mps)/(cfg.c-velocity_mps);
eta_hat = (cfg.c+compensation_velocity_mps)/ ...
    (cfg.c-compensation_velocity_mps);
reference_fft = fft(waveforms);
kernels = complex(zeros(numel(lags),m,m));
peak_matrix = zeros(m,m);
peak_lag_samples = zeros(m,m);
center = floor(n/2)+1;
for i = 1:m
    scaled = interp1(t,waveforms(:,i),eta*t,'linear',0);
    received = scaled.*exp(1j*2*pi*(eta-1)*cfg.fc*t);
    if compensation_velocity_mps ~= 0
        query_time = t/eta_hat;
        compensated = interp1(t,received,query_time,'linear',0) .* ...
            exp(-1j*2*pi*(eta_hat-1)*cfg.fc*query_time);
    else
        compensated = received;
    end
    received_fft = fft(compensated);
    for k = 1:m
        circular = ifft(received_fft.*conj(reference_fft(:,k)));
        shifted = fftshift(circular);
        reference_peak = sum(abs(waveforms(:,k)).^2)+eps;
        cut = shifted((center-max_lag):(center+max_lag))/reference_peak;
        kernels(:,i,k) = cut;
        [peak_matrix(i,k),idx] = max(abs(cut));
        peak_lag_samples(i,k) = lags(idx);
    end
end
end
