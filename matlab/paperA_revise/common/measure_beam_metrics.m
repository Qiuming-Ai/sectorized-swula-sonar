function metrics = measure_beam_metrics(pattern, angles_deg, focus_deg, sector_edges)
%MEASURE_BEAM_METRICS HPBW, sidelobes and sector-edge response.
p = abs(pattern(:).');
p = p / max(p + eps);
p_db = 20*log10(max(p, 1e-12));
[~, peak_idx] = max(p);

left3 = find(p_db(1:peak_idx) <= -3, 1, 'last');
right_rel = find(p_db(peak_idx:end) <= -3, 1, 'first');
if isempty(left3), left3 = 1; end
if isempty(right_rel), right3 = numel(p); else, right3 = peak_idx + right_rel - 1; end

left_min = local_minimum_left(p, peak_idx);
right_min = local_minimum_right(p, peak_idx);
side_mask = true(size(p));
side_mask(left_min:right_min) = false;
if any(side_mask)
    psl_db = 20*log10(max(p(side_mask)) + eps);
else
    psl_db = NaN;
end

edge_gain = interp1(angles_deg, p_db, sector_edges, 'linear', 'extrap');
metrics.focus_deg = focus_deg;
metrics.peak_angle_deg = angles_deg(peak_idx);
metrics.peak_error_deg = metrics.peak_angle_deg - focus_deg;
metrics.hpbw_deg = angles_deg(right3) - angles_deg(left3);
metrics.pslr_db = psl_db;
metrics.left_edge_gain_db = edge_gain(1);
metrics.right_edge_gain_db = edge_gain(2);
metrics.min_edge_gain_db = min(edge_gain);
metrics.peak_gain_linear = max(abs(pattern));
end

function idx = local_minimum_left(p, peak_idx)
idx = 1;
for k = peak_idx-1:-1:2
    if p(k) <= p(k-1) && p(k) <= p(k+1)
        idx = k;
        return;
    end
end
end

function idx = local_minimum_right(p, peak_idx)
idx = numel(p);
for k = peak_idx+1:numel(p)-1
    if p(k) <= p(k-1) && p(k) <= p(k+1)
        idx = k;
        return;
    end
end
end
