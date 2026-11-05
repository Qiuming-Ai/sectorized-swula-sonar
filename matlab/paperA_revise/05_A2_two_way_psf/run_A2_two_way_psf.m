function metrics_table = run_A2_two_way_psf(mode)
%RUN_A2_TWO_WAY_PSF Broadband two-way spatial and range PSF comparison.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A2(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% The stitched H(scan,target)=Tx(target)*Rx(scan,target) exposes grating lobes.
n_f = numel(cfg.frequency_grid);
leakage_score = zeros(n_f,1);
maps = cell(n_f,1);
for jf = 1:n_f
    maps{jf} = stitched_two_way_map(cfg,cfg.frequency_grid(jf));
    leakage_score(jf) = off_sector_score(cfg,maps{jf});
end
[~,worst_idx] = max(leakage_score);
worst_frequency = cfg.frequency_grid(worst_idx);
center_map = stitched_two_way_map(cfg,cfg.fc);
worst_map = maps{worst_idx};

fig1 = figure('Color','w','Position',[80 80 1250 500]);
subplot(1,2,1); imagesc(cfg.angle_grid_deg,cfg.angle_grid_deg, ...
    20*log10(abs(center_map)/(max(abs(center_map(:)))+eps)+eps));
axis xy; clim([-50 0]); colorbar; xlabel('Target angle (deg)');
ylabel('Scan angle (deg)'); title(sprintf('Two-way H at %.0f kHz',cfg.fc/1e3));
subplot(1,2,2); imagesc(cfg.angle_grid_deg,cfg.angle_grid_deg, ...
    20*log10(abs(worst_map)/(max(abs(worst_map(:)))+eps)+eps));
axis xy; clim([-50 0]); colorbar; xlabel('Target angle (deg)');
ylabel('Scan angle (deg)'); title(sprintf('Worst tested frequency: %.0f kHz', ...
    worst_frequency/1e3));
colormap(turbo); sgtitle('Stitched transmit-receive spatial operator');
save_figure_pair(fig1,fullfile(out_dir,'A2_two_way_operator_maps'));

[full_wave,~] = make_waveform_bank(cfg,"fullband");
[full_kernel,~,full_lags] = xcorr_kernel_bank(full_wave,cfg.max_range_lag_samples);
[packed_wave,packed_meta] = make_waveform_bank(cfg,"fdm_current");
[packed_kernel,~,packed_lags] = xcorr_kernel_bank(packed_wave,cfg.max_range_lag_samples);
[fdm_wave,fdm_meta] = make_waveform_bank(cfg,"fdm_strict");
[fdm_kernel,~,fdm_lags] = xcorr_kernel_bank(fdm_wave,cfg.max_range_lag_samples);
[cdm_wave,~] = make_waveform_bank(cfg,"cdm");
[cdm_kernel,~,cdm_lags] = xcorr_kernel_bank(cdm_wave,cfg.max_range_lag_samples);

n_case = numel(cfg.psf_target_angles_deg);
n_method = numel(cfg.psf_methods);
records = repmat(struct(),n_case*n_method,1);
psfs = cell(n_case,n_method);
row = 0;
fig2 = figure('Color','w','Position',[20 40 1750 850]);
for ic = 1:n_case
    target_angle = cfg.psf_target_angles_deg(ic);
    [~,target_idx] = min(abs(cfg.angle_grid_deg-target_angle));
    sector = find_sector(cfg,target_angle);
    for im = 1:n_method
        method = cfg.psf_methods(im);
        if method == "packed_fdm"
            range_cut = abs(packed_kernel(:,sector,sector));
            lags = packed_lags;
        elseif method == "sector_fdm"
            range_cut = abs(fdm_kernel(:,sector,sector));
            lags = fdm_lags;
        elseif method == "cdm_zc"
            range_cut = abs(cdm_kernel(:,sector,sector));
            lags = cdm_lags;
        else
            range_cut = abs(full_kernel(:,1,1));
            lags = full_lags;
        end
        angle_cut = method_angle_psf(cfg,method,target_idx,sector, ...
            packed_meta,fdm_meta);
        psf = range_cut(:)*abs(angle_cut(:)).';
        psfs{ic,im} = psf;
        q = psf_metrics(psf,lags,cfg.angle_grid_deg,target_idx,cfg);
        row = row+1;
        records(row).TargetAngleDeg = target_angle;
        records(row).Method = method;
        records(row).PeakAngleErrorDeg = q.PeakAngleErrorDeg;
        records(row).RangeFWHMM = q.RangeFWHMM;
        records(row).AngleFWHMDeg = q.AngleFWHMDeg;
        records(row).PSLRDb = q.PSLRDb;
        records(row).ISLRDb = q.ISLRDb;
        records(row).GratingPeakDb = q.GratingPeakDb;
        subplot(n_case,n_method,(ic-1)*n_method+im);
        imagesc(cfg.angle_grid_deg,lags*cfg.range_sample_spacing_m, ...
            20*log10(psf/(max(psf(:))+eps)+eps));
        axis xy; clim([-50 0]); xlim([target_angle-20,target_angle+20]);
        xlabel('Angle (deg)'); ylabel('Range offset (m)');
        title(sprintf('%s, target %.1f^o',strrep(method,'_','\_'),target_angle));
    end
end
colormap(turbo); sgtitle('Two-dimensional point-spread functions');
save_figure_pair(fig2,fullfile(out_dir,'A2_method_psf_comparison'));

metrics_table = struct2table(records);
writetable(metrics_table,fullfile(out_dir,'A2_psf_metrics.csv'));
frequency_table = table(cfg.frequency_grid(:)/1e3,leakage_score, ...
    'VariableNames',{'FrequencyKHz','OffSectorScore'});
writetable(frequency_table,fullfile(out_dir,'A2_frequency_leakage.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A2_parameters.json'));
save(fullfile(out_dir,'A2_workspace.mat'),'cfg','metrics_table','frequency_table', ...
    'center_map','worst_map','worst_frequency','psfs');
fprintf('A2 complete. Results: %s\n',out_dir);
disp(metrics_table);
end

function h = stitched_two_way_map(cfg,frequency)
angles = cfg.angle_grid_deg;
[tx,rx] = array_responses(cfg,frequency,angles);
masks = sector_masks(cfg,angles);
h = complex(zeros(numel(angles),numel(angles)));
for k = 1:size(masks,1)
    scan = find(masks(k,:));
    h(scan,:) = rx(scan,:).*tx(k,:);
end
end

function score = off_sector_score(cfg,h)
masks = sector_masks(cfg,cfg.angle_grid_deg);
undesired = [];
desired = [];
for k = 1:size(masks,1)
    desired = [desired;reshape(abs(h(masks(k,:),masks(k,:))).^2,[],1)]; %#ok<AGROW>
    undesired = [undesired;reshape(abs(h(masks(k,:),~masks(k,:))).^2,[],1)]; %#ok<AGROW>
end
score = max(undesired)/(mean(desired(:))+eps);
end

function sector = find_sector(cfg,angle)
sector = find(angle>=cfg.sector_edges_deg(:,1) & ...
    angle<=cfg.sector_edges_deg(:,2),1,'first');
if isempty(sector), error('Target angle is outside the configured FOV.'); end
end

function cut = method_angle_psf(cfg,method,target_idx,sector, ...
    packed_meta,fdm_meta)
angles = cfg.angle_grid_deg;
switch method
    case "dense_lambda2"
        cut = coherent_rx_average(cfg,target_idx,cfg.n_rx_dense, ...
            cfg.dense_rx_spacing_lambda,cfg.frequency_grid);
    case "sw_ula_wideband"
        cut = coherent_rx_average(cfg,target_idx,cfg.n_rx_sparse, ...
            cfg.rx_spacing_lambda,cfg.frequency_grid);
    case "packed_fdm"
        frequency = packed_meta.center_frequency_hz(sector);
        [tx,rx] = array_responses(cfg,frequency,angles);
        cut = tx(sector,target_idx).*rx(:,target_idx).';
        masks = sector_masks(cfg,angles);
        cut(~masks(sector,:)) = 0;
    case "sector_fdm"
        frequency = fdm_meta.center_frequency_hz(sector);
        [tx,rx] = array_responses(cfg,frequency,angles);
        cut = tx(sector,target_idx).*rx(:,target_idx).';
        masks = sector_masks(cfg,angles);
        cut(~masks(sector,:)) = 0;
    case {"cdm_zc","ideal_cdm"}
        cut = complex(zeros(size(angles)));
        for frequency = cfg.frequency_grid
            [tx,rx] = array_responses(cfg,frequency,angles);
            cut = cut + tx(sector,target_idx).*rx(:,target_idx).';
        end
        cut = cut/numel(cfg.frequency_grid);
        masks = sector_masks(cfg,angles);
        cut(~masks(sector,:)) = 0;
    otherwise
        error('Unknown PSF method: %s',method);
end
cut = cut/(max(abs(cut))+eps);
end

function cut = coherent_rx_average(cfg,target_idx,nrx,spacing,frequencies)
cut = complex(zeros(size(cfg.angle_grid_deg)));
for frequency = frequencies
    [~,rx] = array_responses(cfg,frequency,cfg.angle_grid_deg, ...
        Nrx=nrx,RxSpacingLambda=spacing);
    cut = cut + rx(:,target_idx).';
end
cut = cut/numel(frequencies);
end

function q = psf_metrics(psf,lags,angles,target_idx,cfg)
p = abs(psf); p = p/(max(p(:))+eps);
[~,peak_linear] = max(p(:));
[peak_r,peak_a] = ind2sub(size(p),peak_linear);
guard = false(size(p));
[rlo,rhi] = first_minimum_bounds(p(:,peak_a),peak_r,cfg.psf_guard(1));
[alo,ahi] = first_minimum_bounds(p(peak_r,:),peak_a,cfg.psf_guard(2));
rsel = rlo:rhi;
asel = alo:ahi;
guard(rsel,asel)=true;
q.PeakAngleErrorDeg = angles(peak_a)-angles(target_idx);
q.RangeFWHMM = half_power_width(p(:,peak_a),lags)*cfg.range_sample_spacing_m;
q.AngleFWHMDeg = half_power_width_axis(p(peak_r,:),angles);
q.PSLRDb = 20*log10(max(p(~guard))+eps);
q.GratingPeakDb = q.PSLRDb;
q.ISLRDb = 10*log10((sum(p(~guard).^2)+eps)/(sum(p(guard).^2)+eps));
end

function width = half_power_width(x,lags)
x=x(:)/(max(x)+eps); [~,k]=max(x);
l=find(x(1:k)<=1/sqrt(2),1,'last'); r0=find(x(k:end)<=1/sqrt(2),1,'first');
if isempty(l),l=1;end
if isempty(r0),r=numel(x);else,r=k+r0-1;end
width=lags(r)-lags(l);
end

function width = half_power_width_axis(x,axis_values)
x=x(:).'/(max(x)+eps); [~,k]=max(x);
l=find(x(1:k)<=1/sqrt(2),1,'last'); r0=find(x(k:end)<=1/sqrt(2),1,'first');
if isempty(l),l=1;end
if isempty(r0),r=numel(x);else,r=k+r0-1;end
width=axis_values(r)-axis_values(l);
end

function [lo,hi] = first_minimum_bounds(x,peak,fallback)
x = abs(x(:));
lo = max(1,peak-fallback);
hi = min(numel(x),peak+fallback);
for k=peak-1:-1:2
    if x(k)<=x(k-1) && x(k)<=x(k+1), lo=k; break; end
end
for k=peak+1:numel(x)-1
    if x(k)<=x(k-1) && x(k)<=x(k+1), hi=k; break; end
end
end
