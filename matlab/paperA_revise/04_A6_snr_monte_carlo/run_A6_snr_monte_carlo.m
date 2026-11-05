function aggregate = run_A6_snr_monte_carlo(mode)
%RUN_A6_SNR_MONTE_CARLO Paired random-scene/noise Monte-Carlo sweep.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A6(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

models = cell(size(cfg.methods));
for im = 1:numel(cfg.methods)
    models{im} = prepare_reduced_model(cfg,cfg.methods(im));
end

n_total = numel(cfg.scenario_names)*numel(cfg.snr_values_db)* ...
    cfg.monte_carlo_runs*numel(cfg.methods);
records = repmat(struct(),n_total,1);
row = 0;
for is = 1:numel(cfg.scenario_names)
    for irun = 1:cfg.monte_carlo_runs
        scene_seed = cfg.default_seed + is*10000 + irun;
        [scene,~] = generate_scatter_scene(cfg, ...
            cfg.scenario_scatter_counts(is),cfg.scenario_range_layers(is), ...
            cfg.scenario_dynamic_ranges_db(is),scene_seed);
        for isnr = 1:numel(cfg.snr_values_db)
            snr_db = cfg.snr_values_db(isnr);
            for im = 1:numel(cfg.methods)
                row = row+1;
                noise_seed = scene_seed + isnr*100 + im;
                image_data = apply_reduced_model(cfg,scene,models{im},snr_db,noise_seed);
                q = compute_image_metrics(scene,image_data,cfg.metric_guard_radius);
                records(row).Scenario = cfg.scenario_names(is);
                records(row).SNRDb = snr_db;
                records(row).Run = irun;
                records(row).Method = cfg.methods(im);
                records(row).NMSEdB = q.NMSEdB;
                records(row).GlobalSSIM = q.GlobalSSIM;
                records(row).SupportLeakage = q.SupportLeakage;
                records(row).GhostToMainDb = q.GhostToMainDb;
                records(row).CNR = q.CNR;
            end
        end
        fprintf('A6 %s run %d/%d complete.\n',cfg.scenario_names(is), ...
            irun,cfg.monte_carlo_runs);
    end
end
raw = struct2table(records);
writetable(raw,fullfile(out_dir,'A6_all_trials.csv'));
aggregate = aggregate_trials(raw,cfg);
writetable(aggregate,fullfile(out_dir,'A6_summary_mean_ci.csv'));

fig1 = figure('Color','w','Position',[50 50 1500 470]);
for is = 1:numel(cfg.scenario_names)
    subplot(1,numel(cfg.scenario_names),is); hold on;
    for im = 1:numel(cfg.methods)
        rows = aggregate.Scenario==cfg.scenario_names(is) & ...
            aggregate.Method==cfg.methods(im);
        errorbar(aggregate.SNRDb(rows),aggregate.MeanNMSEdB(rows), ...
            aggregate.CI95NMSEdB(rows),'o-','LineWidth',1.2, ...
            'DisplayName',cfg.methods(im));
    end
    grid on; xlabel('Output SNR (dB)'); ylabel('NMSE (dB)');
    title(cfg.scenario_names(is));
    if is == numel(cfg.scenario_names), legend('Location','eastoutside'); end
end
sgtitle(sprintf('Monte-Carlo image error (%d trials per point)', ...
    cfg.monte_carlo_runs));
save_figure_pair(fig1,fullfile(out_dir,'A6_nmse_vs_snr'));

fig2 = figure('Color','w','Position',[70 70 1500 470]);
for is = 1:numel(cfg.scenario_names)
    subplot(1,numel(cfg.scenario_names),is); hold on;
    for im = 1:numel(cfg.methods)
        rows = aggregate.Scenario==cfg.scenario_names(is) & ...
            aggregate.Method==cfg.methods(im);
        errorbar(aggregate.SNRDb(rows),aggregate.MeanCNR(rows), ...
            aggregate.CI95CNR(rows),'o-','LineWidth',1.2, ...
            'DisplayName',cfg.methods(im));
    end
    grid on; xlabel('Output SNR (dB)'); ylabel('CNR');
    title(cfg.scenario_names(is));
end
sgtitle('Monte-Carlo contrast-to-noise ratio');
save_figure_pair(fig2,fullfile(out_dir,'A6_cnr_vs_snr'));

parameter_snapshot(cfg,fullfile(out_dir,'A6_parameters.json'));
save(fullfile(out_dir,'A6_workspace.mat'),'cfg','raw','aggregate');
fprintf('A6 complete. Results: %s\n',out_dir);
end

function aggregate = aggregate_trials(raw,cfg)
rows = repmat(struct(),0,1);
for is = 1:numel(cfg.scenario_names)
    for isnr = 1:numel(cfg.snr_values_db)
        for im = 1:numel(cfg.methods)
            pick = raw.Scenario==cfg.scenario_names(is) & ...
                raw.SNRDb==cfg.snr_values_db(isnr) & raw.Method==cfg.methods(im);
            x_nmse = raw.NMSEdB(pick); x_cnr = raw.CNR(pick);
            x_ssim = raw.GlobalSSIM(pick); x_leak = raw.SupportLeakage(pick);
            n = numel(x_nmse);
            item.Scenario = cfg.scenario_names(is);
            item.SNRDb = cfg.snr_values_db(isnr);
            item.Method = cfg.methods(im);
            item.N = n;
            item.MeanNMSEdB = mean(x_nmse);
            item.CI95NMSEdB = 1.96*std(x_nmse)/sqrt(n);
            item.MeanCNR = mean(x_cnr);
            item.CI95CNR = 1.96*std(x_cnr)/sqrt(n);
            item.MeanGlobalSSIM = mean(x_ssim);
            item.CI95GlobalSSIM = 1.96*std(x_ssim)/sqrt(n);
            item.MeanSupportLeakage = mean(x_leak);
            item.CI95SupportLeakage = 1.96*std(x_leak)/sqrt(n);
            if isempty(rows)
                rows = item;
            else
                rows(end+1,1) = item; %#ok<AGROW>
            end
        end
    end
end
aggregate = struct2table(rows);
end
