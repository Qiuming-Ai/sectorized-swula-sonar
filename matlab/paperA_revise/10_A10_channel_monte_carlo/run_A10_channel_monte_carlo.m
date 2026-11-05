function summary = run_A10_channel_monte_carlo(mode)
%RUN_A10_CHANNEL_MONTE_CARLO Channel-level paired-noise uncertainty study.
if nargin<1, mode="quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A10(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

[scenes,catalogs] = make_scenes(cfg);
models = cell(1,numel(cfg.methods));
for im=1:numel(cfg.methods)
    models{im} = prepare_channel_model(cfg,cfg.methods(im));
end

records = repmat(struct(),numel(cfg.scenario_names)*numel(cfg.snr_values_db)* ...
    cfg.monte_carlo_runs*numel(cfg.methods),1);
row=0;
for is=1:numel(cfg.scenario_names)
    clean = cell(1,numel(cfg.methods));
    forward_diag = cell(1,numel(cfg.methods));
    for im=1:numel(cfg.methods)
        [clean{im},forward_diag{im}] = forward_channel_data(cfg,scenes{is},models{im});
    end
    n_raw = size(clean{1},1);
    for isnr=1:numel(cfg.snr_values_db)
        snr_db = cfg.snr_values_db(isnr);
        for irun=1:cfg.monte_carlo_runs
            seed = cfg.noise_seed_base+is*100000+isnr*1000+irun;
            paired = standardized_complex_noise(n_raw,cfg.n_rx_dense,seed);
            for im=1:numel(cfg.methods)
                row=row+1;
                n_rx=models{im}.n_rx;
                [raw_noisy,noise_diag] = add_paired_channel_noise(clean{im},snr_db,paired(:,1:n_rx));
                image_data = reconstruct_channel_image(cfg,raw_noisy,models{im});
                q = compute_image_metrics(scenes{is},image_data,cfg.metric_guard_radius);
                records(row).Scenario = cfg.scenario_names(is);
                records(row).InputSNRDb = snr_db;
                records(row).Run = irun;
                records(row).NoiseSeed = seed;
                records(row).Method = cfg.methods(im);
                records(row).ActualInputSNRDb = noise_diag.actual_snr_db;
                records(row).RawSignalPower = forward_diag{im}.raw_signal_power;
                records(row).NMSEdB = q.NMSEdB;
                records(row).GlobalSSIM = q.GlobalSSIM;
                records(row).SupportLeakage = q.SupportLeakage;
                records(row).GhostToMainDb = q.GhostToMainDb;
                records(row).CNR = q.CNR;
            end
            fprintf('A10 %s SNR %+g dB: run %d/%d complete.\n', ...
                cfg.scenario_names(is),snr_db,irun,cfg.monte_carlo_runs);
        end
    end
end

raw = struct2table(records);
writetable(raw,fullfile(out_dir,'A10_all_trials.csv'));
summary = aggregate_records(raw,cfg);
writetable(summary,fullfile(out_dir,'A10_summary_mean_ci.csv'));
paired_summary = paired_method_differences(raw,cfg);
writetable(paired_summary,fullfile(out_dir,'A10_paired_differences.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A10_parameters.json'));
save(fullfile(out_dir,'A10_workspace.mat'),'cfg','raw','summary', ...
    'paired_summary','scenes','catalogs','-v7');
fprintf('A10 complete. Results: %s\n',out_dir);
end

function [scenes,catalogs] = make_scenes(cfg)
point_angle = [0,0,0,0,0,0,2,-2,6,-6,12,-12,20,-20,36,-36,42,-42].';
point_range = [4.98,5.00,5.04,5.12,5.28,5.50,5.12,5.12,5.12,5.12, ...
    5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12].';
scene_point = complex(zeros(cfg.n_range_samples,numel(cfg.angle_grid_deg)));
range_idx=zeros(numel(point_range),1); angle_idx=range_idx;
for k=1:numel(point_range)
    [~,range_idx(k)] = min(abs(cfg.physical_range_axis_m-point_range(k)));
    [~,angle_idx(k)] = min(abs(cfg.angle_grid_deg-point_angle(k)));
    scene_point(range_idx(k),angle_idx(k)) = scene_point(range_idx(k),angle_idx(k))+1;
end
catalogs{1}=table(range_idx,angle_idx,point_range,point_angle, ...
    'VariableNames',{'RangeIndex','AngleIndex','RangeM','AngleDeg'});
[scene_dense,catalogs{2}] = generate_scatter_scene(cfg,1000,50,40,cfg.scene_seed);
scenes={scene_point,scene_dense};
end

function summary=aggregate_records(raw,cfg)
items=repmat(struct(),0,1);
metric_names=["NMSEdB","GlobalSSIM","SupportLeakage","GhostToMainDb","CNR"];
for is=1:numel(cfg.scenario_names)
    for isnr=1:numel(cfg.snr_values_db)
        for im=1:numel(cfg.methods)
            pick=raw.Scenario==cfg.scenario_names(is) & ...
                raw.InputSNRDb==cfg.snr_values_db(isnr) & raw.Method==cfg.methods(im);
            item.Scenario=cfg.scenario_names(is);
            item.InputSNRDb=cfg.snr_values_db(isnr);
            item.Method=cfg.methods(im);
            item.N=nnz(pick);
            for q=1:numel(metric_names)
                x=raw.(metric_names(q))(pick);
                item.("Mean"+metric_names(q))=mean(x);
                item.("CI95"+metric_names(q))=1.96*std(x)/sqrt(numel(x));
                item.("SD"+metric_names(q))=std(x);
            end
            if isempty(items)
                items=item;
            else
                items(end+1,1)=item; %#ok<AGROW>
            end
        end
    end
end
summary=struct2table(items);
end

function paired=paired_method_differences(raw,cfg)
items=repmat(struct(),0,1);
comparisons=["cdm","dense";"fdm_current","dense";"cdm","fdm_current"];
for is=1:numel(cfg.scenario_names)
    for isnr=1:numel(cfg.snr_values_db)
        block=raw(raw.Scenario==cfg.scenario_names(is) & ...
            raw.InputSNRDb==cfg.snr_values_db(isnr),:);
        for ic=1:size(comparisons,1)
            a=sortrows(block(block.Method==comparisons(ic,1),:),'Run');
            b=sortrows(block(block.Method==comparisons(ic,2),:),'Run');
            d=a.NMSEdB-b.NMSEdB;
            item.Scenario=cfg.scenario_names(is);
            item.InputSNRDb=cfg.snr_values_db(isnr);
            item.Contrast=comparisons(ic,1)+" minus "+comparisons(ic,2);
            item.N=numel(d);
            item.MeanPairedNMSEDifferenceDb=mean(d);
            item.CI95PairedNMSEDifferenceDb=1.96*std(d)/sqrt(numel(d));
            if isempty(items)
                items=item;
            else
                items(end+1,1)=item; %#ok<AGROW>
            end
        end
    end
end
paired=struct2table(items);
end
