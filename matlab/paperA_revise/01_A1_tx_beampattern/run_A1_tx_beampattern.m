function summary = run_A1_tx_beampattern(mode)
%RUN_A1_TX_BEAMPATTERN Exact 12-element Tx array-factor audit.
if nargin < 1, mode = "quick"; end
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A1(mode);
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

angles = cfg.angle_grid_deg;
freqs = cfg.frequency_grid;
m = numel(cfg.focus_angles_deg);
n_f = numel(freqs);
records = repmat(struct(),n_f*m,1);
coverage_db = zeros(n_f,numel(angles));
adjacent_overlap = zeros(n_f,m-1);
row = 0;

for jf = 1:n_f
    [tx,~,~] = array_responses(cfg,freqs(jf),angles);
    broadside_reference = max(abs(tx(ceil(m/2),:)));
    normalized_beams = abs(tx) ./ max(abs(tx),[],2);
    coverage_db(jf,:) = 20*log10(max(normalized_beams,[],1)+eps);
    for ib = 1:m
        row = row+1;
        bm = measure_beam_metrics(tx(ib,:),angles,cfg.focus_angles_deg(ib), ...
            cfg.sector_edges_deg(ib,:));
        records(row).FrequencyKHz = freqs(jf)/1e3;
        records(row).Beam = ib;
        records(row).FocusDeg = bm.focus_deg;
        records(row).PeakDeg = bm.peak_angle_deg;
        records(row).PeakErrorDeg = bm.peak_error_deg;
        records(row).HPBWDeg = bm.hpbw_deg;
        records(row).PSLRDb = bm.pslr_db;
        records(row).LeftEdgeGainDb = bm.left_edge_gain_db;
        records(row).RightEdgeGainDb = bm.right_edge_gain_db;
        records(row).MinEdgeGainDb = bm.min_edge_gain_db;
        records(row).ScanLossDb = 20*log10((bm.peak_gain_linear+eps) / ...
            (broadside_reference+eps));
    end
    power_beams = normalized_beams.^2;
    for ib = 1:m-1
        numerator = trapz(angles,min(power_beams(ib,:),power_beams(ib+1,:)));
        denominator = min(trapz(angles,power_beams(ib,:)), ...
            trapz(angles,power_beams(ib+1,:)));
        adjacent_overlap(jf,ib) = numerator/(denominator+eps);
    end
end
metrics_table = struct2table(records);
writetable(metrics_table,fullfile(out_dir,'A1_beam_metrics.csv'));

overlap_names = arrayfun(@(a,b)sprintf('Beam%d_%d',a,b), ...
    1:m-1,2:m,'UniformOutput',false);
overlap_table = array2table(adjacent_overlap,'VariableNames',overlap_names);
overlap_table.FrequencyKHz = freqs(:)/1e3;
overlap_table = movevars(overlap_table,'FrequencyKHz','Before',1);
writetable(overlap_table,fullfile(out_dir,'A1_adjacent_overlap.csv'));

[tx_fc,~,~] = array_responses(cfg,cfg.fc,angles);
fig1 = figure('Color','w','Position',[80 80 1050 620]);
hold on;
for ib = 1:m
    p_db = 20*log10(abs(tx_fc(ib,:))/(max(abs(tx_fc(ib,:)))+eps)+eps);
    plot(angles,p_db,'LineWidth',1.5,'DisplayName', ...
        sprintf('Beam %d: %.1f^o',ib,cfg.focus_angles_deg(ib)));
end
for edge = unique(cfg.sector_edges_deg(:)).'
    xline(edge,':','Color',[0.55 0.55 0.55],'HandleVisibility','off');
end
grid on; xlim(cfg.fov_deg+[-5 5]); ylim([-45 2]);
xlabel('Angle (deg)'); ylabel('Normalized magnitude (dB)');
title(sprintf('Exact Tx beampatterns at %.0f kHz',cfg.fc/1e3));
legend('Location','southoutside','NumColumns',3);
save_figure_pair(fig1,fullfile(out_dir,'A1_six_tx_beams_fc'));

fig2 = figure('Color','w','Position',[100 100 1080 520]);
imagesc(angles,freqs/1e3,coverage_db); axis xy;
xlim(cfg.fov_deg); clim([-20 0]); colormap(turbo); colorbar;
xlabel('Angle (deg)'); ylabel('Frequency (kHz)');
title('Best-beam coverage versus angle and frequency (dB)');
save_figure_pair(fig2,fullfile(out_dir,'A1_angle_frequency_coverage'));

fig3 = figure('Color','w','Position',[120 120 900 500]);
plot(freqs/1e3,adjacent_overlap,'o-','LineWidth',1.3);
grid on; xlabel('Frequency (kHz)'); ylabel('Fractional power overlap');
title('Adjacent transmit-beam overlap');
overlap_labels = arrayfun(@(a,b)sprintf('Beam %d-%d',a,b), ...
    1:m-1,2:m,'UniformOutput',false);
legend(overlap_labels,'Location','eastoutside');
save_figure_pair(fig3,fullfile(out_dir,'A1_adjacent_overlap'));

parameter_snapshot(cfg,fullfile(out_dir,'A1_parameters.json'));
summary.WorstSectorEdgeGainDb = min(metrics_table.MinEdgeGainDb);
summary.WorstPeakErrorDeg = max(abs(metrics_table.PeakErrorDeg));
summary.HPBWRangeDeg = [min(metrics_table.HPBWDeg),max(metrics_table.HPBWDeg)];
fov = angles>=min(cfg.sector_edges_deg(:)) & angles<=max(cfg.sector_edges_deg(:));
summary.WorstCoverageDb = min(coverage_db(:,fov),[],'all');
summary.MaxAdjacentOverlap = max(adjacent_overlap(:));
save(fullfile(out_dir,'A1_workspace.mat'),'cfg','metrics_table', ...
    'overlap_table','coverage_db','summary');
fprintf('A1 complete. Results: %s\n',out_dir);
disp(struct2table(summary));
end
