function [kernels, peak_matrix, lags] = xcorr_kernel_bank(waveforms, max_lag)
%XCORR_KERNEL_BANK Matched-filter kernels without Signal Processing Toolbox.
[n,m] = size(waveforms);
max_lag = min(max_lag, n-1);
lags = (-max_lag:max_lag).';
kernels = complex(zeros(numel(lags),m,m));
peak_matrix = zeros(m,m);
for k = 1:m
    auto = conv(waveforms(:,k), conj(flipud(waveforms(:,k))));
    auto_peak = abs(auto(n)) + eps;
    for i = 1:m
        r = conv(waveforms(:,i), conj(flipud(waveforms(:,k))));
        cut = r((n-max_lag):(n+max_lag)) / auto_peak;
        kernels(:,i,k) = cut;
        peak_matrix(i,k) = max(abs(cut));
    end
end
end
