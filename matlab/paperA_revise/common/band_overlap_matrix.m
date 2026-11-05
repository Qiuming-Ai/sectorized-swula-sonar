function overlap = band_overlap_matrix(lower_hz, upper_hz)
%BAND_OVERLAP_MATRIX Pairwise fractional spectral overlap.
m = numel(lower_hz);
overlap = zeros(m,m);
for i = 1:m
    for j = 1:m
        common = max(0, min(upper_hz(i),upper_hz(j)) - ...
            max(lower_hz(i),lower_hz(j)));
        denom = min(upper_hz(i)-lower_hz(i), upper_hz(j)-lower_hz(j));
        overlap(i,j) = common / max(denom, eps);
    end
end
end
