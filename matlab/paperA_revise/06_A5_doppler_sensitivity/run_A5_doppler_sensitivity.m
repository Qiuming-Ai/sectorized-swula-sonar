function metrics_table = run_A5_doppler_sensitivity(mode)
%RUN_A5_DOPPLER_SENSITIVITY Auto/cross ambiguity and image degradation.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A5(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

[scene,~] = generate_scatter_scene(cfg,1000,50,40,cfg.default_seed+55);
n_method = numel(cfg.methods);
n_velocity = numel(cfg.velocities_mps);
records = repmat(struct(),n_method*n_velocity,1);
models = cell(1,n_method);
waveforms = cell(1,n_method);
zero_auto = cell(1,n_method);
row = 0;

for im = 1:n_method
    models{im} = prepare_reduced_model(cfg,cfg.methods(im));
    [waveforms{im},~] = make_waveform_bank(cfg,cfg.methods(im));
    [~,zero_auto{im},~] = doppler_kernel_bank(cfg,waveforms{im},0, ...
        cfg.max_range_lag_samples);
    for iv = 1:n_velocity
        v = cfg.velocities_mps(iv);
        [kernels,peaks,peak_lags] = doppler_kernel_bank(cfg,waveforms{im},v, ...
            cfg.max_range_lag_samples);
        moving_model = models{im};
        moving_model.kernels = kernels;
        moving_model.cross_peak = peaks;
        image_data = apply_reduced_model(cfg,scene,moving_model,Inf,cfg.default_seed+iv);
        q = compute_image_metrics(scene,image_data,cfg.metric_guard_radius);
        desired = diag(peaks);
        reference = diag(zero_auto{im});
        normalized_desired = desired./(reference+eps);
        active = false(size(peaks));
        for k = 1:size(peaks,2)
            active(models{im}.source_sets{k},k) = true;
        end
        active(1:size(active,1)+1:end) = false;
        cross = peaks;
        cross(~active) = NaN;
        lag_diag = diag(peak_lags);
        row = row+1;
        records(row).Method = cfg.methods(im);
        records(row).VelocityMps = v;
        records(row).DopplerScale = (cfg.c+v)/(cfg.c-v);
        records(row).WorstAutoPeakLossDb = 20*log10(min(normalized_desired)+eps);
        records(row).MeanAutoPeakLossDb = 20*log10(mean(normalized_desired)+eps);
        records(row).MaxCrossAmbiguity = max(cross,[],'all','omitnan');
        records(row).MaxCrossAmbiguityDb = 20*log10(records(row).MaxCrossAmbiguity+eps);
        records(row).WorstRangeShiftM = max(abs(lag_diag))*cfg.range_sample_spacing_m;
        records(row).NMSEdB = q.NMSEdB;
        records(row).GlobalSSIM = q.GlobalSSIM;
        records(row).SupportLeakage = q.SupportLeakage;
        records(row).GhostToMainDb = q.GhostToMainDb;
        records(row).CNR = q.CNR;
    end
end
metrics_table = struct2table(records);
writetable(metrics_table,fullfile(out_dir,'A5_doppler_metrics.csv'));

% Correlation-form profiles for a centre-sector waveform.  Curves are
% normalized by that waveform's zero-Doppler autocorrelation peak so that
% Doppler-induced peak loss and range migration are both visible.
profile_records = repmat(struct(),0,1);
profile_methods = ["cdm","fdm_current"];
for ipm = 1:numel(profile_methods)
    im = find(cfg.methods==profile_methods(ipm),1);
    [zero_kernel,~,~] = doppler_kernel_bank(cfg,waveforms{im},0, ...
        cfg.max_range_lag_samples);
    beam_index = cfg.profile_beam;
    reference_peak = max(abs(zero_kernel(:,beam_index,beam_index)))+eps;
    for ipv = 1:numel(cfg.profile_velocities_mps)
        velocity = cfg.profile_velocities_mps(ipv);
        [kernel,~,~] = doppler_kernel_bank(cfg,waveforms{im},velocity, ...
            cfg.max_range_lag_samples);
        curve = abs(kernel(:,beam_index,beam_index))/reference_peak;
        range_offset = (-cfg.max_range_lag_samples:cfg.max_range_lag_samples).' * ...
            cfg.range_sample_spacing_m;
        for ilag = 1:numel(range_offset)
            item.Method = profile_methods(ipm);
            item.VelocityMps = velocity;
            item.Beam = beam_index;
            item.RangeOffsetM = range_offset(ilag);
            item.NormalizedMagnitude = curve(ilag);
            item.MagnitudeDb = 20*log10(curve(ilag)+eps);
            if isempty(profile_records)
                profile_records = item;
            else
                profile_records(end+1,1) = item; %#ok<AGROW>
            end
        end
    end
end
profile_table = struct2table(profile_records);
writetable(profile_table,fullfile(out_dir,'A5_correlation_profiles.csv'));

fig1 = figure('Color','w','Position',[50 50 1450 800]);
quantities = ["WorstAutoPeakLossDb","MaxCrossAmbiguityDb", ...
    "WorstRangeShiftM","NMSEdB"];
ylabels = ["Worst auto peak loss (dB)","Max cross ambiguity (dB)", ...
    "Worst range shift (m)","Image NMSE (dB)"];
for iq = 1:4
    subplot(2,2,iq); hold on;
    for im = 1:n_method
        rows = metrics_table.Method==cfg.methods(im);
        plot(metrics_table.VelocityMps(rows),metrics_table.(quantities(iq))(rows), ...
            'o-','LineWidth',1.4,'DisplayName',cfg.methods(im));
    end
    grid on; xlabel('Radial velocity (m/s)'); ylabel(ylabels(iq));
    if iq==1, legend('Location','best'); end
end
sgtitle('Wideband Doppler sensitivity');
save_figure_pair(fig1,fullfile(out_dir,'A5_doppler_summary'));

% Delay-velocity ambiguity map for the worst CDM cross-code pair.
cdm_idx = find(cfg.methods=="cdm",1);
cross_cube = zeros(2*cfg.max_range_lag_samples+1,n_velocity);
pair = [1,2]; worst = -Inf;
for iv = 1:n_velocity
    [~,peaks,~] = doppler_kernel_bank(cfg,waveforms{cdm_idx}, ...
        cfg.velocities_mps(iv),cfg.max_range_lag_samples);
    tmp=peaks; tmp(1:size(tmp,1)+1:end)=NaN;
    tmp_vector = tmp(:);
    tmp_vector(isnan(tmp_vector)) = -Inf;
    [candidate,linear] = max(tmp_vector);
    if candidate>worst
        worst=candidate; [pair(1),pair(2)]=ind2sub(size(tmp),linear);
    end
end
for iv = 1:n_velocity
    [kernels,~,~] = doppler_kernel_bank(cfg,waveforms{cdm_idx}, ...
        cfg.velocities_mps(iv),cfg.max_range_lag_samples);
    cross_cube(:,iv)=abs(kernels(:,pair(1),pair(2)));
end
fig2 = figure('Color','w','Position',[100 100 900 520]);
imagesc(cfg.velocities_mps,(-cfg.max_range_lag_samples:cfg.max_range_lag_samples)* ...
    cfg.range_sample_spacing_m,20*log10(cross_cube/(max(cross_cube(:))+eps)+eps));
axis xy; clim([-50 0]); colorbar; colormap(turbo);
xlabel('Radial velocity (m/s)'); ylabel('Range offset (m)');
title(sprintf('CDM cross-ambiguity: code %d to matched filter %d',pair(1),pair(2)));
save_figure_pair(fig2,fullfile(out_dir,'A5_cdm_cross_ambiguity'));

parameter_snapshot(cfg,fullfile(out_dir,'A5_parameters.json'));
save(fullfile(out_dir,'A5_workspace.mat'),'cfg','metrics_table','profile_table','scene', ...
    'cross_cube','pair');
fprintf('A5 complete. Results: %s\n',out_dir);
disp(metrics_table);
end
