function summary=run_A13_detection_audit()
%RUN_A13_DETECTION_AUDIT Paired-noise detection and localization audit.
this_dir=fileparts(mfilename('fullpath'));
root=fileparts(this_dir);
addpath(fullfile(root,'common'));
S=load(fullfile(this_dir,'results','A13_workspace.mat'),'cfg','scenes','catalogs');
cfg=S.cfg; scene=S.scenes{1}; truth=S.catalogs{1};
thresholds_db=[-6,-8,-10,-12,-15];
methods=cfg.methods;
models=cell(1,numel(methods)); clean_raw=cell(1,numel(methods));
for im=1:numel(methods)
    models{im}=prepare_channel_model(cfg,methods(im));
    clean_raw{im}=forward_channel_data(cfg,scene,models{im});
end
n_raw=size(clean_raw{1},1);
records=repmat(struct(),0,1);
row=0;
for isnr=1:numel(cfg.snr_values_db)
    snr_db=cfg.snr_values_db(isnr);
    for irun=1:cfg.monte_carlo_runs
        seed=cfg.noise_seed_base+100000+isnr*1000+irun;
        paired=standardized_complex_noise(n_raw,cfg.n_rx_dense,seed);
        for im=1:numel(methods)
            noisy=add_paired_channel_noise(clean_raw{im},snr_db,paired(:,1:models{im}.n_rx));
            image_data=reconstruct_channel_image(cfg,noisy,models{im});
            q=compute_detection_assignment_metrics(image_data,cfg.physical_range_axis_m,cfg.angle_grid_deg,truth,thresholds_db, ...
                'NMSRangeM',0.015,'NMSAngleDeg',1,'GateRangeM',0.04,'GateAngleDeg',2);
            for iq=1:numel(q)
                row=row+1; z=q(iq);
                z.InputSNRDb=snr_db;
                z.Run=irun;
                z.NoiseSeed=seed;
                z.Method=methods(im);
                if row==1, records=z; else, records(row,1)=z; end
            end
        end
        fprintf('A13 detection SNR %+g: %d/%d\n',snr_db,irun,cfg.monte_carlo_runs);
    end
end
raw=struct2table(records(1:row));
summary=aggregate_detection(raw,cfg,thresholds_db);
writetable(raw,fullfile(this_dir,'results','A13_detection_trials.csv'));
writetable(summary,fullfile(this_dir,'results','A13_detection_summary.csv'));
save(fullfile(this_dir,'results','A13_detection_workspace.mat'),'cfg','raw','summary','thresholds_db');
disp(summary);
end

function summary=aggregate_detection(raw,cfg,thresholds_db)
metrics=["DetectionFraction","MeanAngleAbsErrorDeg","RMSEAngleDeg","MeanRangeAbsErrorMm","RMSERangeMm","FalsePeakCount","StrongestUnassignedPeakDb"];
out=repmat(struct(),0,1);
for isnr=1:numel(cfg.snr_values_db)
    for im=1:numel(cfg.methods)
        for it=1:numel(thresholds_db)
            pick=raw.InputSNRDb==cfg.snr_values_db(isnr)&raw.Method==cfg.methods(im)&raw.ThresholdDb==thresholds_db(it);
            z.InputSNRDb=cfg.snr_values_db(isnr);z.Method=cfg.methods(im);z.ThresholdDb=thresholds_db(it);z.N=nnz(pick);
            for k=1:numel(metrics)
                x=raw.(metrics(k))(pick);x=x(isfinite(x));
                z.("Mean"+metrics(k))=mean(x,'omitnan');
                z.("CI95"+metrics(k))=1.96*std(x,'omitnan')/sqrt(numel(x));
            end
            if isempty(out),out=z;else,out(end+1,1)=z;end %#ok<AGROW>
        end
    end
end
summary=struct2table(out);
end
