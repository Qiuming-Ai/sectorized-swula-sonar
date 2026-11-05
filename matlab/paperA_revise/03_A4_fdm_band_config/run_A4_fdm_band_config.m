function summary_table = run_A4_fdm_band_config(mode)
%RUN_A4_FDM_BAND_CONFIG Compare current, strict and optimized FDM layouts.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A4(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

if cfg.optimize_assignment
    [permutation,cost_table] = optimize_slot_assignment(cfg);
    cfg.fdm_optimized_permutation = permutation;
    writetable(cost_table,fullfile(out_dir,'A4_optimized_assignment_costs.csv'));
else
    permutation = cfg.fdm_optimized_permutation;
end
fprintf('Optimized strict-slot permutation: [%s]\n',num2str(permutation));

n_method = numel(cfg.methods);
models = cell(1,n_method);
waveforms = cell(1,n_method);
wave_meta = cell(1,n_method);
overlaps = cell(1,n_method);
cross_peaks = cell(1,n_method);
raw_cross_peaks = cell(1,n_method);
records = repmat(struct(),n_method,1);
psf = cell(1,n_method);
psf_lags = cell(1,n_method);

for im = 1:n_method
    [waveforms{im},wave_meta{im}] = make_waveform_bank(cfg,cfg.methods(im));
    [~,raw_cross_peaks{im}] = xcorr_kernel_bank(waveforms{im}, ...
        size(waveforms{im},1)-1);
    [kernels,cross_peaks{im},psf_lags{im}] = ...
        xcorr_kernel_bank(waveforms{im},cfg.max_range_lag_samples);
    overlaps{im} = band_overlap_matrix(wave_meta{im}.lower_edge_hz, ...
        wave_meta{im}.upper_edge_hz);
    demultiplexing = "matched filtering";
    models{im} = prepare_reduced_model(cfg,cfg.methods(im));
    desired_psfs = zeros(size(kernels,1),size(kernels,2));
    for k = 1:size(kernels,2)
        desired_psfs(:,k) = abs(kernels(:,k,k));
    end
    psf{im} = mean(desired_psfs,2);
    desired_psf = psf{im};
    records(im).Method = cfg.methods(im);
    records(im).OccupiedLowerKHz = min(wave_meta{im}.lower_edge_hz)/1e3;
    records(im).OccupiedUpperKHz = max(wave_meta{im}.upper_edge_hz)/1e3;
    records(im).MaxFractionalOverlap = active_off_diagonal_max(overlaps{im},models{im});
    records(im).MaxProcessedCrossCoupling = ...
        active_off_diagonal_max(cross_peaks{im},models{im});
    records(im).MaxRawLinearCrossCorrelation = ...
        active_off_diagonal_max(raw_cross_peaks{im},models{im});
    records(im).RangeFWHMSamples = fwhm_samples(desired_psf);
    records(im).RangeFWHMM = records(im).RangeFWHMSamples*cfg.range_sample_spacing_m;
    records(im).WorstCouplingDb = model_coupling_db(models{im});
    records(im).Demultiplexing = demultiplexing;
end
summary_table = struct2table(records);
writetable(summary_table,fullfile(out_dir,'A4_fdm_configuration_metrics.csv'));

fig1 = figure('Color','w','Position',[50 50 1300 820]);
for im = 1:n_method
    subplot(n_method,1,im); hold on;
    [n,m] = size(waveforms{im});
    nfft = 2^nextpow2(4*n);
    f = ((-nfft/2:nfft/2-1)/nfft)*cfg.waveform_fs/1e3 + cfg.fc/1e3;
    for k = 1:m
        spectrum = fftshift(abs(fft(waveforms{im}(:,k),nfft)));
        plot(f,20*log10(spectrum/(max(spectrum)+eps)+eps),'LineWidth',1.0);
    end
    xlim([390 510]); ylim([-60 3]); grid on;
    ylabel('Magnitude (dB)'); title(strrep(cfg.methods(im),'_','\_'));
end
xlabel('Frequency (kHz)');
sgtitle('FDM waveform spectra');
save_figure_pair(fig1,fullfile(out_dir,'A4_fdm_spectra'));

fig2 = figure('Color','w','Position',[80 80 1200 420]);
for im = 1:n_method
    subplot(1,n_method,im);
    display_overlap = nan(size(overlaps{im}));
    for k = 1:size(display_overlap,2)
        display_overlap(models{im}.source_sets{k},k) = ...
            overlaps{im}(models{im}.source_sets{k},k);
    end
    imagesc(display_overlap,[0 1]); axis image; colorbar;
    xlabel('Waveform'); ylabel('Waveform');
    title(strrep(cfg.methods(im),'_','\_'));
end
sgtitle('Simultaneous-channel fractional passband overlap');
save_figure_pair(fig2,fullfile(out_dir,'A4_overlap_matrices'));

fig3 = figure('Color','w','Position',[100 100 900 520]); hold on;
for im = 1:n_method
    plot(psf_lags{im}*cfg.range_sample_spacing_m, ...
        20*log10(psf{im}/max(psf{im})+eps), ...
        'LineWidth',1.5,'DisplayName',cfg.methods(im));
end
grid on; ylim([-60 2]); xlabel('Range offset (m)'); ylabel('Mean MF response (dB)');
legend('Location','best'); title('Range point-spread functions');
save_figure_pair(fig3,fullfile(out_dir,'A4_range_psf'));

[scene,~] = generate_scatter_scene(cfg,200,10,20,cfg.default_seed+44);
images = cell(1,n_method);
fig4 = figure('Color','w','Position',[40 40 1400 480]);
subplot(1,n_method+1,1); plot_image_db(cfg,scene,50,'Truth');
for im = 1:n_method
    images{im} = apply_reduced_model(cfg,scene,models{im},Inf,cfg.default_seed+44);
    subplot(1,n_method+1,im+1);
    plot_image_db(cfg,images{im},50,cfg.methods(im));
end
sgtitle('Paired medium-complexity scene');
save_figure_pair(fig4,fullfile(out_dir,'A4_scene_comparison'));

parameter_snapshot(cfg,fullfile(out_dir,'A4_parameters.json'));
save(fullfile(out_dir,'A4_workspace.mat'),'cfg','summary_table','overlaps', ...
    'cross_peaks','raw_cross_peaks','wave_meta','psf','psf_lags', ...
    'scene','images','permutation');
fprintf('A4 complete. Results: %s\n',out_dir);
disp(summary_table);
end

function [best_perm,cost_table] = optimize_slot_assignment(cfg)
m = numel(cfg.focus_angles_deg);
sub_bw = (1-cfg.fdm_guard_fraction)*cfg.total_bw/m;
slot_width = cfg.total_bw/m;
slots = (-cfg.total_bw/2+slot_width/2)+(0:m-1)*slot_width;
angles = cfg.angle_grid_deg;
masks = sector_masks(cfg,angles);
cost = zeros(m,m);
for ib = 1:m
    own = masks(ib,:);
    for is = 1:m
        [tx,rx] = array_responses(cfg,cfg.fc+slots(is),angles);
        h = abs(rx(own,:).*tx(ib,:)).^2;
        desired = mean(h(:,own),'all')+eps;
        cost(ib,is) = mean(h(:,~own),'all')/desired;
    end
end
all_perm = perms(1:m);
scores = zeros(size(all_perm,1),1);
for ip = 1:size(all_perm,1)
    idx = sub2ind([m,m],1:m,all_perm(ip,:));
    scores(ip) = sum(cost(idx));
end
[~,best] = min(scores);
best_perm = all_perm(best,:);
cost_table = array2table(cost,'VariableNames',compose('Slot%d',1:m));
cost_table.Beam = (1:m).';
cost_table = movevars(cost_table,'Beam','Before',1);
cost_table.AssignedSlot = best_perm(:);
cost_table.SubbandBandwidthKHz = repmat(sub_bw/1e3,m,1);
end

function width = fwhm_samples(x)
x = abs(x(:))/max(abs(x(:))+eps);
[~,peak] = max(x);
left = find(x(1:peak)<=0.5,1,'last');
right_rel = find(x(peak:end)<=0.5,1,'first');
if isempty(left), left=1; end
if isempty(right_rel), right=numel(x); else, right=peak+right_rel-1; end
width = right-left;
end

function db = model_coupling_db(model)
m = size(model.spatial,1);
ratios = [];
for k = 1:m
    desired = norm(model.spatial{k,k},'fro')^2+eps;
    for i = model.source_sets{k}
        if i==k, continue; end
        spatial_ratio = norm(model.spatial{k,i},'fro')^2/desired;
        ratios(end+1) = spatial_ratio*model.cross_peak(i,k)^2; %#ok<AGROW>
    end
end
db = 10*log10(max(ratios)+eps);
end

function value = active_off_diagonal_max(matrix,model)
values = [];
for k = 1:size(matrix,2)
    source_set = model.source_sets{k};
    source_set(source_set==k) = [];
    values = [values;reshape(matrix(source_set,k),[],1)]; %#ok<AGROW>
end
if isempty(values), value=0; else, value=max(values); end
end
