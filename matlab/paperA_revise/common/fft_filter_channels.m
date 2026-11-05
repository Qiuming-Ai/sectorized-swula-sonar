function y = fft_filter_channels(x, h, n_output)
%FFT_FILTER_CHANNELS Linear convolution along the first matrix dimension.

if nargin < 3
    n_output = size(x,1)+numel(h)-1;
end
n_linear = size(x,1)+numel(h)-1;
if n_output > n_linear
    error('n_output cannot exceed the linear-convolution length.');
end
n_fft = 2^nextpow2(n_linear);
H = fft(h(:),n_fft);
Y = ifft(fft(x,n_fft,1).*H,[],1);
y = Y(1:n_output,:);
end
