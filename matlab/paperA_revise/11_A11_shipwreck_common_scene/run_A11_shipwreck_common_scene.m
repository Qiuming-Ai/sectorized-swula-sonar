function metrics = run_A11_shipwreck_common_scene(mode,geometry_dir,arrival_source,out_dir,mex_dir)
%RUN_A11_SHIPWRECK_COMMON_SCENE Fair STL-driven extended-scatter comparison.
if nargin<1, mode="quick"; end
if nargin<2, geometry_dir=[]; end
if nargin<3, arrival_source="archive"; end
this_dir=fileparts(mfilename('fullpath'));
root=fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg=config_A11(mode,geometry_dir,arrival_source);
if nargin<4 || isempty(out_dir), out_dir=cfg.output_dir; end
if nargin>=5 && ~isempty(mex_dir), cfg.cuda_mex_dir=mex_dir; end
cfg.output_dir=out_dir;
if ~exist(out_dir,'dir'), mkdir(out_dir); end

[scene,visual_reference,ship_catalog,seabed_catalog,arrival_catalog,geometry]= ...
    build_shipwreck_scene(cfg);
n_method=numel(cfg.methods);
models=cell(1,n_method); clean_raw=cell(1,n_method);
images_clean=complex(zeros(cfg.n_range_samples,numel(cfg.angle_grid_deg),n_method));
images_noisy=complex(zeros(size(images_clean)));
forward_diag=cell(1,n_method); noise_diag=cell(1,n_method);
for im=1:n_method
    models{im}=prepare_channel_model(cfg,cfg.methods(im));
    [clean_raw{im},forward_diag{im}]=forward_channel_data(cfg,scene,models{im});
end

n_raw=size(clean_raw{1},1);
paired=standardized_complex_noise(n_raw,cfg.n_rx_dense,cfg.noise_seed);
for im=1:n_method
    images_clean(:,:,im)=reconstruct_channel_image(cfg,clean_raw{im},models{im});
    [raw_noisy,noise_diag{im}]=add_paired_channel_noise(clean_raw{im}, ...
        cfg.input_snr_db,paired(:,1:models{im}.n_rx));
    images_noisy(:,:,im)=reconstruct_channel_image(cfg,raw_noisy,models{im});
    fprintf('A11 %s reconstruction complete.\n',cfg.methods(im));
end

reference=images_noisy(:,:,1);
items=repmat(struct(),n_method,1);
for im=1:n_method
    qt=compute_image_metrics(abs(scene),images_noisy(:,:,im),cfg.metric_guard_radius);
    qr=compute_image_metrics(abs(reference),images_noisy(:,:,im),cfg.metric_guard_radius);
    ref_range=sum(abs(reference),2); est_range=sum(abs(images_noisy(:,:,im)),2);
    ref_angle=sum(abs(reference),1).'; est_angle=sum(abs(images_noisy(:,:,im)),1).';
    items(im).Method=cfg.methods(im);
    items(im).InputSNRDb=noise_diag{im}.actual_snr_db;
    items(im).NMSEToTruthDb=qt.NMSEdB;
    items(im).GSSIMToTruth=qt.GlobalSSIM;
    items(im).CNRToTruthSupport=qt.CNR;
    items(im).TruthSupportLeakage=qt.SupportLeakage;
    items(im).NMSEToDenseDb=qr.NMSEdB;
    items(im).GSSIMToDense=qr.GlobalSSIM;
    items(im).RangeProfileCorrelation=corr_scalar(ref_range,est_range);
    items(im).AngleProfileCorrelation=corr_scalar(ref_angle,est_angle);
    items(im).RawSignalPower=forward_diag{im}.raw_signal_power;
end
metrics=struct2table(items);
writetable(metrics,fullfile(out_dir,'A11_shipwreck_metrics.csv'));
writetable(ship_catalog,fullfile(out_dir,'A11_ship_facets.csv'));
writetable(seabed_catalog,fullfile(out_dir,'A11_seabed_scatterers.csv'));
writetable(arrival_catalog,fullfile(out_dir,'A11_raytrace_arrivals.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A11_parameters.json'));

angle_deg=cfg.angle_grid_deg;
range_m=cfg.physical_range_axis_m;
truth_magnitude=visual_reference;
image_magnitude=abs(images_noisy);
method_labels=cfg.method_labels;
save(fullfile(out_dir,'A11_plot_data.mat'),'angle_deg','range_m', ...
    'truth_magnitude','image_magnitude','method_labels','geometry','-v7');
save(fullfile(out_dir,'A11_workspace.mat'),'cfg','scene','visual_reference', ...
    'ship_catalog','seabed_catalog','arrival_catalog','geometry','images_clean','images_noisy', ...
    'metrics','forward_diag','noise_diag','-v7.3');
fprintf('A11 complete. Results: %s\n',out_dir);
disp(metrics);
end

function value=corr_scalar(a,b)
a=a(:); b=b(:);
a=a-mean(a); b=b-mean(b);
value=sum(a.*b)/sqrt((sum(a.^2)+eps)*(sum(b.^2)+eps));
end
