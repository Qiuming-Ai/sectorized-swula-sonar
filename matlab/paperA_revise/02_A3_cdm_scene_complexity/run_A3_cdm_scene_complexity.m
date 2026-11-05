function results = run_A3_cdm_scene_complexity(mode)
%RUN_A3_CDM_SCENE_COMPLEXITY Quantify failure as scene support grows.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A3(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

conditions = make_conditions(cfg);
models = cell(size(cfg.methods));
for im = 1:numel(cfg.methods)
    fprintf('Preparing %s model...\n',cfg.methods(im));
    models{im} = prepare_reduced_model(cfg,cfg.methods(im));
end

n_row = height(conditions)*numel(cfg.methods);
records = repmat(struct(),n_row,1);
row = 0;
for ic = 1:height(conditions)
    seed = cfg.default_seed + ic;
    [scene,~] = generate_scatter_scene(cfg,conditions.ScatterCount(ic), ...
        conditions.RangeLayers(ic),conditions.DynamicRangeDb(ic),seed);
    for im = 1:numel(cfg.methods)
        row = row+1;
        image_data = apply_reduced_model(cfg,scene,models{im},Inf,seed);
        q = compute_image_metrics(scene,image_data,cfg.metric_guard_radius);
        records(row).Condition = string(conditions.Condition(ic));
        records(row).Factor = string(conditions.Factor(ic));
        records(row).FactorValue = conditions.FactorValue(ic);
        records(row).ScatterCount = conditions.ScatterCount(ic);
        records(row).RangeLayers = conditions.RangeLayers(ic);
        records(row).DynamicRangeDb = conditions.DynamicRangeDb(ic);
        records(row).Method = cfg.methods(im);
        records(row).NMSEdB = q.NMSEdB;
        records(row).GlobalSSIM = q.GlobalSSIM;
        records(row).SupportLeakage = q.SupportLeakage;
        records(row).GhostToMainDb = q.GhostToMainDb;
        records(row).CNR = q.CNR;
    end
    fprintf('A3 condition %d/%d complete.\n',ic,height(conditions));
end
results = struct2table(records);
writetable(results,fullfile(out_dir,'A3_scene_complexity_metrics.csv'));

fig1 = figure('Color','w','Position',[60 60 1250 900]);
factors = ["ScatterCount","RangeLayers","DynamicRangeDb"];
for jf = 1:numel(factors)
    subplot(3,1,jf); hold on;
    subset = results(results.Factor==factors(jf),:);
    for im = 1:numel(cfg.methods)
        rows = subset.Method==cfg.methods(im);
        plot(subset.FactorValue(rows),subset.NMSEdB(rows),'o-', ...
            'LineWidth',1.4,'DisplayName',cfg.methods(im));
    end
    grid on; xlabel(factors(jf)); ylabel('NMSE (dB)');
    if jf == 1, legend('Location','eastoutside'); end
end
sgtitle('Image degradation versus scene complexity');
save_figure_pair(fig1,fullfile(out_dir,'A3_complexity_nmse'));

% Save a common difficult scene so the visual comparison is exactly paired.
[hard_scene,hard_catalog] = generate_scatter_scene(cfg,max(cfg.scatter_counts), ...
    max(cfg.range_layer_counts),max(cfg.dynamic_ranges_db),cfg.default_seed+999);
hard_images = cell(size(cfg.methods));
hard_records = repmat(struct(),numel(cfg.methods),1);
fig2 = figure('Color','w','Position',[30 50 1750 480]);
subplot(1,numel(cfg.methods)+1,1);
plot_image_db(cfg,hard_scene,50,'Truth');
for im = 1:numel(cfg.methods)
    hard_images{im} = apply_reduced_model(cfg,hard_scene,models{im},Inf, ...
        cfg.default_seed+999);
    q = compute_image_metrics(hard_scene,hard_images{im},cfg.metric_guard_radius);
    hard_records(im).Method = cfg.methods(im);
    hard_records(im).NMSEdB = q.NMSEdB;
    hard_records(im).GlobalSSIM = q.GlobalSSIM;
    hard_records(im).SupportLeakage = q.SupportLeakage;
    hard_records(im).GhostToMainDb = q.GhostToMainDb;
    hard_records(im).CNR = q.CNR;
    subplot(1,numel(cfg.methods)+1,im+1);
    plot_image_db(cfg,hard_images{im},50,cfg.methods(im));
end
sgtitle('Paired 1000-scatter / 50-layer / 40-dB scene');
save_figure_pair(fig2,fullfile(out_dir,'A3_hard_scene_comparison'));
hard_metrics = struct2table(hard_records);
writetable(hard_metrics,fullfile(out_dir,'A3_hard_1000scatter_50layer_metrics.csv'));

parameter_snapshot(cfg,fullfile(out_dir,'A3_parameters.json'));
save(fullfile(out_dir,'A3_workspace.mat'),'cfg','conditions','results', ...
    'hard_scene','hard_catalog','hard_images','hard_metrics');
fprintf('A3 complete. Results: %s\n',out_dir);
end

function conditions = make_conditions(cfg)
if cfg.full_factorial
    [s,l,d] = ndgrid(cfg.scatter_counts,cfg.range_layer_counts,cfg.dynamic_ranges_db);
    s=s(:); l=l(:); d=d(:);
    factor = repmat("FullFactorial",numel(s),1);
    value = (1:numel(s)).';
else
    s = cfg.scatter_counts(:);
    l = min(s,repmat(cfg.reference_layer_count,numel(s),1));
    d = repmat(cfg.reference_dynamic_range_db,numel(s),1);
    factor = repmat("ScatterCount",numel(s),1); value = s;

    l2 = cfg.range_layer_counts(:);
    s2 = repmat(cfg.reference_scatter_count,numel(l2),1);
    d2 = repmat(cfg.reference_dynamic_range_db,numel(l2),1);
    s=[s;s2]; l=[l;l2]; d=[d;d2];
    factor=[factor;repmat("RangeLayers",numel(l2),1)]; value=[value;l2];

    d3 = cfg.dynamic_ranges_db(:);
    s3 = repmat(cfg.reference_scatter_count,numel(d3),1);
    l3 = repmat(cfg.reference_layer_count,numel(d3),1);
    s=[s;s3]; l=[l;l3]; d=[d;d3];
    factor=[factor;repmat("DynamicRangeDb",numel(d3),1)]; value=[value;d3];
end
condition = factor + "_" + string(value);
conditions = table(condition,factor,value,s,l,d,'VariableNames', ...
    {'Condition','Factor','FactorValue','ScatterCount','RangeLayers','DynamicRangeDb'});
end
