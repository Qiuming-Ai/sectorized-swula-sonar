function [waveforms, meta] = make_waveform_bank(cfg, family)
%MAKE_WAVEFORM_BANK Generate original-paper LFM and ZC baseband waveforms.
family = lower(string(family));
fs = cfg.waveform_fs;
duration_s = cfg.pulse_duration_s;
n = max(32, round(duration_s * fs));
t = (0:n-1).' / fs;
tau = t - mean(t);
win = local_window("hamming", n);
m = numel(cfg.focus_angles_deg);

switch family
    case {"fullband", "dense"}
        centers = zeros(1,m);
        bandwidths = cfg.total_bw * ones(1,m);
    case "cdm"
        centers = zeros(1,m);
        bandwidths = cfg.total_bw * ones(1,m);
    case "fdm_strict"
        sub_bw = (1-cfg.fdm_guard_fraction) * cfg.total_bw/m;
        slot_width = cfg.total_bw/m;
        centers = (-cfg.total_bw/2 + slot_width/2) + (0:m-1)*slot_width;
        bandwidths = sub_bw * ones(1,m);
    case "fdm_current"
        centers = cfg.current_fdm_centers_hz - cfg.fc;
        bandwidths = cfg.current_fdm_bandwidth_hz * ones(1,m);
    case "fdm_optimized"
        sub_bw = (1-cfg.fdm_guard_fraction) * cfg.total_bw/m;
        slot_width = cfg.total_bw/m;
        slots = (-cfg.total_bw/2 + slot_width/2) + (0:m-1)*slot_width;
        centers = slots(cfg.fdm_optimized_permutation);
        bandwidths = sub_bw * ones(1,m);
    otherwise
        error('Unknown waveform family: %s', family);
end

waveforms = complex(zeros(n,m));
for k = 1:m
    bw = bandwidths(k);
    if family == "cdm"
        l = cfg.zc_length;
        chip_n = (0:l-1).';
        root = cfg.zc_roots(k);
        if gcd(root,l) ~= 1
            error('ZC root %d is not coprime with code length %d.',root,l);
        end
        zc = exp(-1j*pi*root.*chip_n.*(chip_n+1)/l);
        samples_per_chip = round(fs/cfg.zc_chip_rate_hz);
        if abs(samples_per_chip-fs/cfg.zc_chip_rate_hz)>1e-12
            error('waveform_fs / zc_chip_rate_hz must be an integer.');
        end
        s = repelem(zc,samples_per_chip);
        if numel(s)<n
            s(end+1:n,1) = 0;
        else
            s = s(1:n);
        end
    else
        chirp = exp(1j*pi*(bw/duration_s).*tau.^2 + ...
            1j*2*pi*centers(k).*tau);
        s = win .* chirp;
    end
    waveforms(:,k) = s / sqrt(sum(abs(s).^2) + eps);
end

meta.family = family;
meta.fs_hz = fs;
meta.time_s = t;
meta.center_offset_hz = centers;
meta.center_frequency_hz = cfg.fc + centers;
meta.bandwidth_hz = bandwidths;
meta.lower_edge_hz = meta.center_frequency_hz - bandwidths/2;
meta.upper_edge_hz = meta.center_frequency_hz + bandwidths/2;
meta.duration_s = duration_s;
if family == "cdm"
    meta.code_family = "Zadoff-Chu";
    meta.code_length = cfg.zc_length;
    meta.code_roots = cfg.zc_roots;
    meta.chip_rate_hz = cfg.zc_chip_rate_hz;
end
meta.slot_index = ones(1,m);
meta.n_slots = max(meta.slot_index);
end
