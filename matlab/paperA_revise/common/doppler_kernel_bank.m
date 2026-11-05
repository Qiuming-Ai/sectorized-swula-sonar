function [kernels,peak_matrix,peak_lag_samples] = doppler_kernel_bank(cfg,waveforms,velocity_mps,max_lag)
%DOPPLER_KERNEL_BANK Wideband time-scale Doppler matched-filter kernels.
% Positive velocity denotes approach. The analytic baseband signal is time
% scaled and receives the residual carrier shift after demodulation.
[n,m] = size(waveforms);
max_lag = min(max_lag,n-1);
t = (0:n-1).'/cfg.waveform_fs;
eta = (cfg.c+velocity_mps)/(cfg.c-velocity_mps);
lags = (-max_lag:max_lag).';
kernels = complex(zeros(numel(lags),m,m));
peak_matrix = zeros(m,m);
peak_lag_samples = zeros(m,m);
for i = 1:m
    scaled = interp1(t,waveforms(:,i),eta*t,'linear',0);
    received = scaled .* exp(1j*2*pi*(eta-1)*cfg.fc*t);
    for k = 1:m
        reference_auto = conv(waveforms(:,k),conj(flipud(waveforms(:,k))));
        reference_peak = abs(reference_auto(n))+eps;
        r = conv(received,conj(flipud(waveforms(:,k))));
        cut = r((n-max_lag):(n+max_lag))/reference_peak;
        kernels(:,i,k) = cut;
        [peak_matrix(i,k),idx] = max(abs(cut));
        peak_lag_samples(i,k) = lags(idx);
    end
end
end
