function summary=run_A16_shipwreck_phase_ensemble(mode,geometry_dir,out_dir,saved_a11_workspace)
%RUN_A16_SHIPWRECK_PHASE_ENSEMBLE Reconstruct 30 paired scattering phases.
if nargin<1,mode="quick";end
if nargin<2,geometry_dir=[];end
this_dir=fileparts(mfilename('fullpath'));root=fileparts(this_dir);addpath(fullfile(root,'common'));
cfg=config_A16(mode,geometry_dir);
if nargin<3 || isempty(out_dir),out_dir=cfg.output_dir;end
cfg.output_dir=out_dir;
if nargin>=4 && ~isempty(saved_a11_workspace) && strlength(string(saved_a11_workspace))>0
    cfg.saved_a11_workspace=saved_a11_workspace;
    saved=load(cfg.saved_a11_workspace,'arrival_catalog','visual_reference');
    arrival=saved.arrival_catalog;visual_reference=saved.visual_reference;
else
    [~,visual_reference,~,~,arrival]=build_shipwreck_scene(cfg);
end
support=false(cfg.n_range_samples,numel(cfg.angle_grid_deg));
if ~isequal(size(visual_reference),size(support))
    error('A16:GridMismatch','Use the full paper-mode A11 archive grid, not the quick grid.');
end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
linear=sub2ind(size(support),arrival.RangeIndex,arrival.AngleIndex);support(linear)=true;
models=cell(1,numel(cfg.methods));for im=1:numel(cfg.methods),models{im}=prepare_channel_model(cfg,cfg.methods(im));end
records=repmat(struct(),cfg.phase_runs*numel(cfg.methods),1);sensitivity=repmat(struct(),cfg.phase_runs*numel(cfg.methods)*size(cfg.guard_radii,1),1);
row=0;srow=0;representative_images=[];
for irun=1:cfg.phase_runs
    phase_seed=cfg.phase_seed_base+irun-1;rng(phase_seed,'twister');phase=2*pi*rand(height(arrival),1);
    coefficient=arrival.MetricAmplitude.*exp(1j*phase);
    scene=accumarray([arrival.RangeIndex,arrival.AngleIndex],coefficient,size(support),@sum,complex(0));scene=scene/(max(abs(scene(:)))+eps);
    clean_raw=cell(1,numel(cfg.methods));images=complex(zeros(cfg.n_range_samples,numel(cfg.angle_grid_deg),numel(cfg.methods)));
    for im=1:numel(cfg.methods),clean_raw{im}=forward_channel_data(cfg,scene,models{im});end
    n_raw=size(clean_raw{1},1);noise_seed=cfg.noise_seed_base+irun;paired=standardized_complex_noise(n_raw,cfg.n_rx_dense,noise_seed);
    noise_diag=cell(1,numel(cfg.methods));
    for im=1:numel(cfg.methods)
        [raw_noisy,noise_diag{im}]=add_paired_channel_noise(clean_raw{im},cfg.input_snr_db,paired(:,1:models{im}.n_rx));
        images(:,:,im)=reconstruct_channel_image(cfg,raw_noisy,models{im});
    end
    dense=images(:,:,1);dense_range=sum(abs(dense),2);dense_angle=sum(abs(dense),1).';
    for im=1:numel(cfg.methods)
        q=compute_image_metrics_v2(visual_reference,images(:,:,im),support,[],CoreRadius=[4,1],GuardRadius=cfg.primary_guard_radius);
        estimate_range=sum(abs(images(:,:,im)),2);estimate_angle=sum(abs(images(:,:,im)),1).';
        row=row+1;records(row).Run=irun;records(row).PhaseSeed=phase_seed;records(row).NoiseSeed=noise_seed;records(row).Method=cfg.methods(im);
        records(row).ActualInputSNRDb=noise_diag{im}.actual_snr_db;records(row).Correlation2DToDense=correlation2d(abs(dense),abs(images(:,:,im)));
        records(row).RangeProfileCorrelationToDense=correlation1d(dense_range,estimate_range);records(row).AngleProfileCorrelationToDense=correlation1d(dense_angle,estimate_angle);
        records(row).GSSIMToTruthEnvelope=q.GlobalSSIM;records(row).Correlation2DToTruthEnvelope=q.TruthCorrelation2D;records(row).PrimarySupportLeakage=q.SupportLeakage;
        for ig=1:size(cfg.guard_radii,1)
            guard=cfg.guard_radii(ig,:);qs=compute_image_metrics_v2(visual_reference,images(:,:,im),support,[],CoreRadius=[4,1],GuardRadius=guard);
            srow=srow+1;sensitivity(srow).Run=irun;sensitivity(srow).Method=cfg.methods(im);sensitivity(srow).GuardRangeBins=guard(1);sensitivity(srow).GuardAngleBins=guard(2);sensitivity(srow).SupportLeakage=qs.SupportLeakage;sensitivity(srow).MaxOutsideDb=qs.MaxOutsideDb;
        end
    end
    if irun==1,representative_images=images;end
    fprintf('A16 phase realization %d/%d complete.\n',irun,cfg.phase_runs);
end
raw=struct2table(records(1:row));sensitivity_table=struct2table(sensitivity(1:srow));summary=aggregate(raw,cfg);
writetable(raw,fullfile(out_dir,'A16_all_trials.csv'));writetable(summary,fullfile(out_dir,'A16_summary_mean_ci.csv'));writetable(sensitivity_table,fullfile(out_dir,'A16_leakage_sensitivity.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A16_parameters.json'));save(fullfile(out_dir,'A16_workspace.mat'),'cfg','raw','summary','sensitivity_table','representative_images','support','-v7.3');disp(summary);
end
function value=correlation1d(a,b),a=a(:)-mean(a(:));b=b(:)-mean(b(:));value=sum(a.*b)/sqrt((sum(a.^2)+eps)*(sum(b.^2)+eps));end
function value=correlation2d(a,b),value=correlation1d(a(:),b(:));end
function summary=aggregate(raw,cfg)
names=["Correlation2DToDense","RangeProfileCorrelationToDense","AngleProfileCorrelationToDense","GSSIMToTruthEnvelope","Correlation2DToTruthEnvelope","PrimarySupportLeakage"];
items=repmat(struct(),numel(cfg.methods),1);
for im=1:numel(cfg.methods)
 pick=raw.Method==cfg.methods(im);items(im).Method=cfg.methods(im);items(im).N=nnz(pick);
 for iq=1:numel(names),x=raw.(names(iq))(pick);items(im).("Mean"+names(iq))=mean(x);items(im).("CI95"+names(iq))=1.96*std(x)/sqrt(numel(x));end
end
summary=struct2table(items);
end
