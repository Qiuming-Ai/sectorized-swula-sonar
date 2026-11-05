function summary = run_A12_fdm_ambiguity_audit()
%RUN_A12_FDM_AMBIGUITY_AUDIT Compare full-delay and processing-window peaks.
this_dir = fileparts(mfilename('fullpath'));
root = fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg = config_A12();
out_dir = fullfile(this_dir,'results');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

pair_items = repmat(struct(),0,1);
summary_items = repmat(struct(),0,1);
profile_blocks = cell(numel(cfg.methods)*numel(cfg.velocities_mps),1);
spectrum_blocks = cell(numel(cfg.methods)*numel(cfg.focus_angles_deg),1);
profile_block_index = 0;
spectrum_block_index = 0;
worst_pair_by_method = zeros(numel(cfg.methods),2);

for im = 1:numel(cfg.methods)
    method = cfg.methods(im);
    [waveforms,meta] = make_waveform_bank(cfg,method);
    n = size(waveforms,1);
    full_half = n-1;
    lags = (-full_half:full_half).';
    delay_ms = 1e3*lags/cfg.waveform_fs;
    range_offset_m = cfg.c*lags/(2*cfg.waveform_fs);

    % Identify the worst zero-Doppler ordered pair on the complete pulse.
    [zero_kernels,~,~] = doppler_kernel_bank(cfg,waveforms,0,full_half);
    zero_peak = -Inf;
    worst_pair = [1,2];
    for i = 1:size(waveforms,2)
        for k = 1:size(waveforms,2)
            if i==k, continue; end
            candidate = max(abs(zero_kernels(:,i,k)));
            if candidate > zero_peak
                zero_peak = candidate;
                worst_pair = [i,k];
            end
        end
    end
    worst_pair_by_method(im,:) = worst_pair;

    for iv = 1:numel(cfg.velocities_mps)
        velocity = cfg.velocities_mps(iv);
        [kernels,~,~] = doppler_kernel_bank(cfg,waveforms,velocity,full_half);
        local_mask = abs(lags)<=cfg.processing_half_window_samples;
        full_worst = -Inf;
        local_worst = -Inf;
        zero_worst = -Inf;
        full_pair = [NaN,NaN];
        local_pair = [NaN,NaN];

        for i = 1:size(waveforms,2)
            for k = 1:size(waveforms,2)
                if i==k, continue; end
                curve = abs(kernels(:,i,k));
                [full_value,full_idx] = max(curve);
                local_indices = find(local_mask);
                [local_value,local_rel_idx] = max(curve(local_mask));
                local_idx = local_indices(local_rel_idx);
                zero_idx = full_half+1;
                zero_value = curve(zero_idx);

                item.Method = method;
                item.VelocityMps = velocity;
                item.SourceWaveform = i;
                item.MatchedFilter = k;
                item.SourceSectorAngleDeg = cfg.focus_angles_deg(i);
                item.OutputSectorAngleDeg = cfg.focus_angles_deg(k);
                item.SourceCenterKHz = meta.center_frequency_hz(i)/1e3;
                item.OutputCenterKHz = meta.center_frequency_hz(k)/1e3;
                item.AbsoluteSpacingKHz = abs(meta.center_frequency_hz(i)-meta.center_frequency_hz(k))/1e3;
                item.FullPeakDb = 20*log10(full_value+eps);
                item.FullPeakDelayMs = delay_ms(full_idx);
                item.FullPeakRangeOffsetM = range_offset_m(full_idx);
                item.ProcessingWindowPeakDb = 20*log10(local_value+eps);
                item.ProcessingWindowPeakDelayMs = delay_ms(local_idx);
                item.ProcessingWindowPeakRangeOffsetM = range_offset_m(local_idx);
                item.ZeroDelayDb = 20*log10(zero_value+eps);
                pair_items = append_item(pair_items,item);

                if full_value>full_worst
                    full_worst=full_value; full_pair=[i,k];
                end
                if local_value>local_worst
                    local_worst=local_value; local_pair=[i,k];
                end
                zero_worst=max(zero_worst,zero_value);
            end
        end

        s.Method = method;
        s.VelocityMps = velocity;
        s.WorstFullPeakDb = 20*log10(full_worst+eps);
        s.WorstFullSource = full_pair(1);
        s.WorstFullFilter = full_pair(2);
        s.WorstProcessingWindowPeakDb = 20*log10(local_worst+eps);
        s.WorstProcessingWindowSource = local_pair(1);
        s.WorstProcessingWindowFilter = local_pair(2);
        s.WorstZeroDelayDb = 20*log10(zero_worst+eps);
        summary_items = append_item(summary_items,s);

        % Export a complete-delay profile for the method's zero-Doppler worst pair.
        pair = worst_pair_by_method(im,:);
        curve = abs(kernels(:,pair(1),pair(2)));
        profile_block_index = profile_block_index+1;
        profile_blocks{profile_block_index} = table( ...
            repmat(method,numel(lags),1), ...
            repmat(velocity,numel(lags),1), ...
            repmat(pair(1),numel(lags),1), ...
            repmat(pair(2),numel(lags),1), ...
            lags,delay_ms,range_offset_m,curve,20*log10(curve+eps),local_mask, ...
            'VariableNames',{'Method','VelocityMps','SourceWaveform', ...
            'MatchedFilter','LagSamples','DelayMs','RangeOffsetM','Magnitude', ...
            'MagnitudeDb','InProcessingWindow'});
    end

    % Export normalized spectra using the same finite waveforms.
    f = (-cfg.spectrum_fft_length/2:cfg.spectrum_fft_length/2-1).' * ...
        cfg.waveform_fs/cfg.spectrum_fft_length;
    for k = 1:size(waveforms,2)
        spectrum = abs(fftshift(fft(waveforms(:,k),cfg.spectrum_fft_length)));
        spectrum = spectrum/(max(spectrum)+eps);
        keep = abs(f)<=0.65*cfg.total_bw;
        selected_f = f(keep);
        selected_spectrum = spectrum(keep);
        spectrum_block_index = spectrum_block_index+1;
        spectrum_blocks{spectrum_block_index} = table( ...
            repmat(method,numel(selected_f),1), ...
            repmat(k,numel(selected_f),1), ...
            repmat(cfg.focus_angles_deg(k),numel(selected_f),1), ...
            selected_f/1e3,(selected_f+cfg.fc)/1e3, ...
            max(-100,20*log10(selected_spectrum+eps)), ...
            'VariableNames',{'Method','Waveform','SectorAngleDeg', ...
            'BasebandFrequencyKHz','RFFrequencyKHz','MagnitudeDb'});
    end
end

pair_table = struct2table(pair_items);
summary = struct2table(summary_items);
profile_table = vertcat(profile_blocks{1:profile_block_index});
spectrum_table = vertcat(spectrum_blocks{1:spectrum_block_index});
writetable(pair_table,fullfile(out_dir,'A12_pairwise_peaks.csv'));
writetable(summary,fullfile(out_dir,'A12_doppler_summary.csv'));
writetable(profile_table,fullfile(out_dir,'A12_delay_profiles.csv'));
writetable(spectrum_table,fullfile(out_dir,'A12_spectra.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A12_parameters.json'));
save(fullfile(out_dir,'A12_workspace.mat'),'cfg','pair_table','summary', ...
    'profile_table','spectrum_table','worst_pair_by_method');
fprintf('A12 complete. Results: %s\n',out_dir);
disp(summary);
end

function records = append_item(records,item)
if isempty(records)
    records = item;
else
    records(end+1,1) = item; %#ok<AGROW>
end
end
