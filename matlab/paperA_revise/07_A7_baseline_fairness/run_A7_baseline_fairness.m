function outputs = run_A7_baseline_fairness(mode)
%RUN_A7_BASELINE_FAIRNESS Same-aperture DAS plus conditional MVDR/MUSIC.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A7(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

geometries = make_geometries(cfg);
names = fieldnames(geometries);
angles = cfg.angle_grid_deg;

geometry_records = repmat(struct(),numel(names),1);
point_rows = repmat(struct(),numel(names)*numel(cfg.point_target_angles_deg),1);
point_row = 0;
broadband_operators = cell(size(names));
for ig = 1:numel(names)
    p = geometries.(names{ig});
    broadband_operators{ig} = broadband_das_operator(p,cfg.frequency_grid,cfg.fc,angles);
    spacing = diff(sort(p));
    geometry_records(ig).Geometry = string(names{ig});
    geometry_records(ig).ElementCount = numel(p);
    geometry_records(ig).ApertureLambda = max(p)-min(p);
    geometry_records(ig).MinimumSpacingLambda = min(spacing);
    geometry_records(ig).MedianSpacingLambda = median(spacing);
    geometry_records(ig).IsDenseOracle = string(names{ig})=="dense_oracle";
    for it = 1:numel(cfg.point_target_angles_deg)
        [~,target_idx] = min(abs(angles-cfg.point_target_angles_deg(it)));
        cut = abs(broadband_operators{ig}(:,target_idx)).';
        cut = cut/(max(cut)+eps);
        [~,peak] = max(cut);
        guard = abs(angles-cfg.point_target_angles_deg(it))<=2;
        point_row = point_row+1;
        point_rows(point_row).Geometry = string(names{ig});
        point_rows(point_row).TargetAngleDeg = cfg.point_target_angles_deg(it);
        point_rows(point_row).PeakAngleDeg = angles(peak);
        point_rows(point_row).LocalizationErrorDeg = angles(peak)- ...
            cfg.point_target_angles_deg(it);
        point_rows(point_row).GhostPeakDb = 20*log10(max(cut(~guard))+eps);
    end
end
geometry_table = struct2table(geometry_records);
point_table = struct2table(point_rows);
writetable(geometry_table,fullfile(out_dir,'A7_geometry_summary.csv'));
writetable(point_table,fullfile(out_dir,'A7_point_das_metrics.csv'));

fig1 = figure('Color','w','Position',[70 70 1000 600]); hold on;
[~,zero_idx] = min(abs(angles));
for ig = 1:numel(names)
    cut = abs(broadband_operators{ig}(:,zero_idx));
    cut = cut/(max(cut)+eps);
    plot(angles,20*log10(cut+eps),'LineWidth',1.3,'DisplayName',names{ig});
end
grid on; xlim([-60 60]); ylim([-55 2]); xlabel('Scan angle (deg)');
ylabel('DAS point response (dB)'); title('Same-aperture broadband DAS, broadside target');
legend('Location','southoutside','NumColumns',3,'Interpreter','none');
save_figure_pair(fig1,fullfile(out_dir,'A7_broadband_das_point_response'));

% Adaptive methods are evaluated only in their valid finite-source setting.
adaptive_rows = repmat(struct(),numel(names)*2,1);
adaptive_spectra = cell(numel(names),2);
row = 0;
fig2 = figure('Color','w','Position',[40 40 1450 800]);
for ig = 1:numel(names)
    p = geometries.(names{ig});
    [mvdr,music] = adaptive_spectra_for_geometry(p,cfg);
    adaptive_spectra{ig,1}=mvdr; adaptive_spectra{ig,2}=music;
    algorithms = ["MVDR","MUSIC"];
    spectra = {mvdr,music};
    for ia = 1:2
        estimates = separated_peaks(angles,spectra{ia}, ...
            numel(cfg.adaptive_source_angles_deg),3);
        row=row+1;
        adaptive_rows(row).Geometry=string(names{ig});
        adaptive_rows(row).Algorithm=algorithms(ia);
        adaptive_rows(row).MeanLocalizationErrorDeg = assignment_error( ...
            cfg.adaptive_source_angles_deg,estimates);
        adaptive_rows(row).EstimatedAnglesDeg = join(string(round(estimates,3)),',');
    end
    subplot(2,ceil(numel(names)/2),ig); hold on;
    plot(angles,10*log10(mvdr/max(mvdr)+eps),'LineWidth',1.2,'DisplayName','MVDR');
    plot(angles,10*log10(music/max(music)+eps),'LineWidth',1.2,'DisplayName','MUSIC');
    for source = cfg.adaptive_source_angles_deg, xline(source,':k','HandleVisibility','off'); end
    ylim([-50 2]); xlim([-60 60]); grid on; title(names{ig},'Interpreter','none');
    xlabel('Angle (deg)'); ylabel('Spectrum (dB)');
    if ig==1, legend('Location','best'); end
end
sgtitle('Three known point sources: MVDR/MUSIC conditional comparison');
save_figure_pair(fig2,fullfile(out_dir,'A7_mvdr_music_point_sources'));
adaptive_table = struct2table(adaptive_rows);
writetable(adaptive_table,fullfile(out_dir,'A7_adaptive_point_metrics.csv'));

% Diffuse-scene comparison uses DAS only; MUSIC source-count assumptions fail.
rng(cfg.default_seed+77,'twister');
x = complex(zeros(numel(angles),1));
idx = randi(numel(angles),cfg.diffuse_scatter_count,1);
amp = 10.^(-30*rand(cfg.diffuse_scatter_count,1)/20).* ...
    exp(1j*2*pi*rand(cfg.diffuse_scatter_count,1));
for k=1:numel(idx), x(idx(k))=x(idx(k))+amp(k); end
dense_idx = find(strcmp(names,'dense_oracle'),1);
dense_image = abs(broadband_operators{dense_idx}*x);
dense_image = dense_image/(max(dense_image)+eps);
diffuse_rows = repmat(struct(),numel(names),1);
fig3 = figure('Color','w','Position',[80 80 1100 600]); hold on;
for ig=1:numel(names)
    estimate = abs(broadband_operators{ig}*x);
    estimate = estimate/(max(estimate)+eps);
    diffuse_rows(ig).Geometry=string(names{ig});
    diffuse_rows(ig).NMSEToDense = sum((estimate-dense_image).^2)/(sum(dense_image.^2)+eps);
    diffuse_rows(ig).CorrelationToDense = corr_toolbox_free(estimate,dense_image);
    plot(angles,20*log10(estimate+eps),'LineWidth',1.1,'DisplayName',names{ig});
end
grid on; ylim([-55 2]); xlabel('Angle (deg)'); ylabel('Normalized DAS image (dB)');
title('100-scatter diffuse scene: DAS only');
legend('Location','southoutside','NumColumns',3,'Interpreter','none');
save_figure_pair(fig3,fullfile(out_dir,'A7_diffuse_scene_das'));
diffuse_table = struct2table(diffuse_rows);
writetable(diffuse_table,fullfile(out_dir,'A7_diffuse_das_metrics.csv'));

parameter_snapshot(cfg,fullfile(out_dir,'A7_parameters.json'));
save(fullfile(out_dir,'A7_workspace.mat'),'cfg','geometries','geometry_table', ...
    'point_table','adaptive_table','diffuse_table','adaptive_spectra');
outputs.Geometry=geometry_table; outputs.PointDAS=point_table;
outputs.Adaptive=adaptive_table; outputs.DiffuseDAS=diffuse_table;
fprintf('A7 complete. Results: %s\n',out_dir);
end

function geometries = make_geometries(cfg)
aperture = (cfg.n_rx_sparse-1)*cfg.rx_spacing_lambda;
geometries.dense_oracle = linspace(-aperture/2,aperture/2,cfg.n_rx_dense).';
geometries.sw_ula = linspace(-aperture/2,aperture/2,cfg.n_rx_sparse).';
rng(cfg.default_seed,'twister');
grid = linspace(-aperture/2,aperture/2,cfg.n_rx_dense);
middle = randperm(cfg.n_rx_dense-2,cfg.n_rx_sparse-2)+1;
geometries.random = sort([grid(1),grid(middle),grid(end)]).';
n1=16; n2=16;
nested = unique([0:n1-1,n1*(1:n2)]);
geometries.nested = rescale_aperture(nested,aperture);
m=16; n=17;
coprime = unique([m*(0:n-1),n*(0:m-1)]);
geometries.coprime = rescale_aperture(coprime,aperture);
end

function p = rescale_aperture(indices,aperture)
p = indices(:)-mean([min(indices),max(indices)]);
p = p/(max(p)-min(p))*aperture;
end

function b = broadband_das_operator(positions_lambda,frequencies,fc,angles)
b = complex(zeros(numel(angles),numel(angles)));
for f=frequencies
    a=exp(1j*2*pi*(f/fc)*positions_lambda(:)*sind(angles));
    b=b+(a'*a)/numel(positions_lambda);
end
b=b/numel(frequencies);
end

function [mvdr,music] = adaptive_spectra_for_geometry(p,cfg)
a_scan=exp(1j*2*pi*p(:)*sind(cfg.angle_grid_deg));
a_src=exp(1j*2*pi*p(:)*sind(cfg.adaptive_source_angles_deg));
signal_cov=a_src*diag(cfg.adaptive_source_powers)*a_src';
signal_power=real(trace(signal_cov))/numel(p);
noise_power=signal_power/10^(cfg.adaptive_snr_db/10);
r=signal_cov+noise_power*eye(numel(p));
loaded=r+1e-8*trace(r)/numel(p)*eye(numel(p));
rinv_a=loaded\a_scan;
mvdr=real(1./(sum(conj(a_scan).*rinv_a,1)+eps));
[v,d]=eig((r+r')/2,'vector');
[~,order]=sort(real(d),'descend'); v=v(:,order);
k=numel(cfg.adaptive_source_angles_deg);
en=v(:,k+1:end);
music=real(1./(sum(abs(en'*a_scan).^2,1)+eps));
mvdr=mvdr(:).'; music=music(:).';
end

function estimates = separated_peaks(angles,spectrum,n_peak,min_separation_deg)
s=spectrum(:).'; estimates=zeros(1,n_peak); available=true(size(s));
for k=1:n_peak
    candidate=s; candidate(~available)=-Inf;
    [~,idx]=max(candidate); estimates(k)=angles(idx);
    available(abs(angles-angles(idx))<min_separation_deg)=false;
end
estimates=sort(estimates);
end

function error_value = assignment_error(truth,estimate)
p=perms(1:numel(truth)); costs=zeros(size(p,1),1);
for k=1:size(p,1), costs(k)=mean(abs(truth-estimate(p(k,:)))); end
error_value=min(costs);
end

function r = corr_toolbox_free(x,y)
x=x(:)-mean(x(:)); y=y(:)-mean(y(:));
r=real(sum(conj(x).*y)/(sqrt(sum(abs(x).^2)*sum(abs(y).^2))+eps));
end
