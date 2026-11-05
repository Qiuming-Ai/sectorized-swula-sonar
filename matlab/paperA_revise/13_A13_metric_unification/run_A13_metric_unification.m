function summary=run_A13_metric_unification(mode)
%RUN_A13_METRIC_UNIFICATION One metric implementation for all main results.
if nargin<1, mode="quick"; end
this_dir=fileparts(mfilename('fullpath'));
root=fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg=config_A13(mode);
out_dir=fullfile(this_dir,'results');
if ~exist(out_dir,'dir'),mkdir(out_dir);end
[scenes,catalogs]=make_scenes(cfg);
models=cell(1,numel(cfg.methods));
for im=1:numel(cfg.methods),models{im}=prepare_channel_model(cfg,cfg.methods(im));end

n_record=numel(cfg.scenario_names)*numel(cfg.snr_values_db)*cfg.monte_carlo_runs*numel(cfg.methods);
records=repmat(struct(),n_record,1); row=0;
clean_images=cell(numel(cfg.scenario_names),numel(cfg.methods));
for is=1:numel(cfg.scenario_names)
    support=support_from_catalog(catalogs{is},cfg.n_range_samples,numel(cfg.angle_grid_deg));
    clean_raw=cell(1,numel(cfg.methods));
    forward_diag=cell(1,numel(cfg.methods));
    for im=1:numel(cfg.methods)
        [clean_raw{im},forward_diag{im}]=forward_channel_data(cfg,scenes{is},models{im});
        clean_images{is,im}=reconstruct_channel_image(cfg,clean_raw{im},models{im});
    end
    n_raw=size(clean_raw{1},1);
    for isnr=1:numel(cfg.snr_values_db)
        snr_db=cfg.snr_values_db(isnr);
        for irun=1:cfg.monte_carlo_runs
            seed=cfg.noise_seed_base+is*100000+isnr*1000+irun;
            paired=standardized_complex_noise(n_raw,cfg.n_rx_dense,seed);
            for im=1:numel(cfg.methods)
                [raw_noisy,noise_diag]=add_paired_channel_noise(clean_raw{im},snr_db,paired(:,1:models{im}.n_rx));
                image_data=reconstruct_channel_image(cfg,raw_noisy,models{im});
                q=compute_image_metrics_v2(scenes{is},image_data,support,clean_images{is,im}, ...
                    CoreRadius=cfg.metric_core_radius,GuardRadius=cfg.metric_guard_radius);
                row=row+1;
                records(row).Scenario=cfg.scenario_names(is);
                records(row).InputSNRDb=snr_db;
                records(row).Run=irun;
                records(row).NoiseSeed=seed;
                records(row).Method=cfg.methods(im);
                records(row).ActualInputSNRDb=noise_diag.actual_snr_db;
                records(row).RawSignalPower=forward_diag{im}.raw_signal_power;
                records(row).TruthNMSEdB=q.TruthNMSEdB;
                records(row).NoiseNMSEdB=q.NoiseNMSEdB;
                records(row).TruthCorrelation2D=q.TruthCorrelation2D;
                records(row).Correlation2DToDenseReference=correlation2d( ...
                    abs(clean_images{is,1}),abs(image_data));
                records(row).GlobalSSIM=q.GlobalSSIM;
                records(row).CNR=q.CNR;
                records(row).SupportLeakage=q.SupportLeakage;
                records(row).MaxOutsideDb=q.MaxOutsideDb;
                records(row).CorePixelCount=q.CorePixelCount;
                records(row).GuardPixelCount=q.GuardPixelCount;
                records(row).BackgroundPixelCount=q.BackgroundPixelCount;
            end
            fprintf('A13 %s SNR %+g: %d/%d\n',cfg.scenario_names(is),snr_db,irun,cfg.monte_carlo_runs);
        end
    end
end
raw=struct2table(records(1:row));
summary=aggregate(raw,cfg);
writetable(raw,fullfile(out_dir,'A13_all_trials.csv'));
writetable(summary,fullfile(out_dir,'A13_summary_mean_ci.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A13_parameters.json'));
save(fullfile(out_dir,'A13_workspace.mat'),'cfg','raw','summary','scenes','catalogs','clean_images','-v7.3');
disp(summary);
end

function [scenes,catalogs]=make_scenes(cfg)
point_angle=[0,0,0,0,0,0,2,-2,6,-6,12,-12,20,-20,36,-36,42,-42].';
point_range=[4.98,5.00,5.04,5.12,5.28,5.50,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12].';
scene_point=complex(zeros(cfg.n_range_samples,numel(cfg.angle_grid_deg)));
range_idx=zeros(numel(point_range),1);angle_idx=range_idx;
for k=1:numel(point_range)
    [~,range_idx(k)]=min(abs(cfg.physical_range_axis_m-point_range(k)));
    [~,angle_idx(k)]=min(abs(cfg.angle_grid_deg-point_angle(k)));
    scene_point(range_idx(k),angle_idx(k))=scene_point(range_idx(k),angle_idx(k))+1;
end
catalogs{1}=table(range_idx,angle_idx,point_range,point_angle,'VariableNames',{'RangeIndex','AngleIndex','RangeM','AngleDeg'});
[scene_dense,catalogs{2}]=generate_scatter_scene(cfg,1000,50,40,cfg.scene_seed);
scenes={scene_point,scene_dense};
end

function summary=aggregate(raw,cfg)
metric_names=["TruthNMSEdB","NoiseNMSEdB","TruthCorrelation2D", ...
    "Correlation2DToDenseReference","GlobalSSIM","CNR", ...
    "SupportLeakage","MaxOutsideDb"];
items=repmat(struct(),0,1);
for is=1:numel(cfg.scenario_names)
    for isnr=1:numel(cfg.snr_values_db)
        for im=1:numel(cfg.methods)
            pick=raw.Scenario==cfg.scenario_names(is)&raw.InputSNRDb==cfg.snr_values_db(isnr)&raw.Method==cfg.methods(im);
            item.Scenario=cfg.scenario_names(is);item.InputSNRDb=cfg.snr_values_db(isnr);item.Method=cfg.methods(im);item.N=nnz(pick);
            for iq=1:numel(metric_names)
                x=raw.(metric_names(iq))(pick);x=x(isfinite(x));
                item.("Mean"+metric_names(iq))=mean(x);
                item.("CI95"+metric_names(iq))=1.96*std(x)/sqrt(numel(x));
            end
            if isempty(items),items=item;else,items(end+1,1)=item;end %#ok<AGROW>
        end
    end
end
summary=struct2table(items);
end

function value=correlation2d(a,b)
a=a(:)-mean(a(:));
b=b(:)-mean(b(:));
value=sum(a.*b)/sqrt((sum(a.^2)+eps)*(sum(b.^2)+eps));
end
