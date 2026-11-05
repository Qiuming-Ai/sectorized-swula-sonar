function noise = standardized_complex_noise(n_row,n_col,seed)
%STANDARDIZED_COMPLEX_NOISE Deterministic unit-power circular complex noise.
rng(seed,'twister');
noise = (randn(n_row,n_col)+1j*randn(n_row,n_col))/sqrt(2);
noise = noise/sqrt(mean(abs(noise(:)).^2)+eps);
end
