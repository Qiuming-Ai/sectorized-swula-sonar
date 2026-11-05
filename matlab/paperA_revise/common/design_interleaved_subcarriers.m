function [allocation,diagnostics] = design_interleaved_subcarriers(cfg)
%DESIGN_INTERLEAVED_SUBCARRIERS Deterministic six-way disjoint tone allocation.
% Every beam receives 80 of 480 orthogonal 200-Hz bins spread across 96 kHz.
% The random search minimizes the worst autocorrelation sidelobe over the
% configured range window; a persistent cache makes repeated calls cheap.
persistent cached_allocation cached_diagnostics cached_signature
signature = [cfg.interleaved_n_tones_total,cfg.interleaved_n_beams, ...
    cfg.interleaved_symbol_duration_s,cfg.waveform_fs,cfg.default_seed];
if ~isempty(cached_signature) && isequal(signature,cached_signature)
    allocation = cached_allocation;
    diagnostics = cached_diagnostics;
    return;
end

m = cfg.interleaved_n_beams;
n_total = cfg.interleaved_n_tones_total;
if mod(n_total,m)~=0 || mod(n_total,2)~=0
    error('Total subcarrier count must be even and divisible by beam count.');
end
n_per_beam = n_total/m;
half = n_total/2;
bins = [-half:-1,1:half];
n_symbol = round(cfg.interleaved_symbol_duration_s*cfg.waveform_fs);
lags = (-min(cfg.interleaved_design_max_lag,floor(n_symbol/2)-1): ...
         min(cfg.interleaved_design_max_lag,floor(n_symbol/2)-1)).';
phase_table = exp(1j*2*pi*(lags/n_symbol)*bins);
main_guard = abs(lags)<=cfg.interleaved_psf_guard_samples;

state = rng;
rng(cfg.default_seed+811,'twister');
best_score = Inf;
best = [];
for trial = 1:cfg.interleaved_allocation_trials
    shuffled = bins(randperm(n_total));
    candidate = reshape(shuffled,n_per_beam,m).';
    spans = max(candidate,[],2)-min(candidate,[],2);
    centers = abs(mean(candidate,2));
    if min(spans)<0.90*n_total || max(centers)>0.05*n_total
        continue;
    end
    worst = 0;
    mean_side_power = 0;
    for k = 1:m
        [~,idx] = ismember(candidate(k,:),bins);
        response = abs(sum(phase_table(:,idx),2))/n_per_beam;
        side = response(~main_guard);
        worst = max(worst,max(side));
        mean_side_power = mean_side_power+mean(side.^2)/m;
    end
    score = worst+0.15*sqrt(mean_side_power)+0.002*max(centers);
    if score<best_score
        best_score = score;
        best = sort(candidate,2);
    end
end
rng(state);
if isempty(best)
    error('No feasible interleaved subcarrier allocation was found.');
end

allocation = best;
diagnostics.Score = best_score;
diagnostics.BinSpacingHz = 1/cfg.interleaved_symbol_duration_s;
diagnostics.SpanHz = (max(best,[],2)-min(best,[],2))*diagnostics.BinSpacingHz;
diagnostics.CenterOffsetHz = mean(best,2)*diagnostics.BinSpacingHz;
diagnostics.PSLRDb = zeros(m,1);
for k = 1:m
    [~,idx] = ismember(best(k,:),bins);
    response = abs(sum(phase_table(:,idx),2))/n_per_beam;
    diagnostics.PSLRDb(k) = 20*log10(max(response(~main_guard))+eps);
end
cached_signature = signature;
cached_allocation = allocation;
cached_diagnostics = diagnostics;
end
