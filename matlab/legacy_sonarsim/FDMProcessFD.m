function timing_data = FDMProcessFD(input_file, max_pings, show_plots)
% Legacy desktop FDM workload. Optional historical plot follows the timer.
legacy_sonarsim_setup();
if nargin < 1 || isempty(input_file)
    paths = legacy_desktop_timing_inputs();
    input_file = paths.fdm;
end
if nargin < 2, max_pings = Inf; end
if nargin < 3, show_plots = false; end
[data, count] = legacy_load_timing_fixture(input_file, "FDM", max_pings);
    win = hamming(data.attributes.receive_array_num);
    cfg.rx_xyz = data.attributes.receive_array_position;
    cfg.c = data.attributes.sound_velocity;
    cfg.t0 = data.attributes.sample_delay-data.attributes.pulse_duration/2;                      % 用你的实际 t_start
    cfg.fs = data.attributes.sampling_frequency/data.attributes.decimate_factor/data.attributes.second_decimate_factor;
    % cfg.angles = data.attributes.scan_angle * pi/180;   % 扫描角（rad）
    cfg.fd.nfft   =  2^nextpow2(round(length(data.pings{1}(:,1))/data.attributes.second_decimate_factor));          % 或 64 / 256，越大分数时延越精细，但更慢
    cfg.fd.sign   = 1;

    %% Temp Para
    fc = data.attributes.center_frequency;
    f0 =  mean(data.attributes.center_frequency);
    fs_decimate = (data.attributes.sampling_frequency/data.attributes.decimate_factor);
    % 初始化时间数组
    timing_data = zeros(count, 1);
    for i = 1:count
    % for i = 1:1
        ping_data = data.pings{i};
        tic
        for j = 1:data.attributes.sector_num
        %% Match Filter
        sig_mf = zeros(size(ping_data));
        for m = 1:data.attributes.receive_array_num
            tmp = conv(ping_data(:,m), data.attributes.match_filter_data(:,j), 'same');
            sig_mf(:,m) = tmp*win(m);
        end
        %% Mix & Decimate
        sig_mf_mix = sig_mf.*repmat(exp(1j*2*pi*(fc(j)-f0)*[0:length(sig_mf(:,1))-1]/fs_decimate)',[1 data.attributes.receive_array_num]);
        clear sig_mf;
        sig_mf_fir = filter(data.attributes.second_fir,1,sig_mf_mix);
        clear sig_mf_mix;
        sig_mf_dec = resample(sig_mf_fir,1,data.attributes.second_decimate_factor);
        clear sig_mf_fir;
        %% TVG
        % To-do
        index = data.attributes.scan_angle>=data.attributes.sector_div(j)&data.attributes.scan_angle<data.attributes.sector_div(j+1);
        cfg.angles =  data.attributes.scan_angle(index)/180*pi;
        cfg.BW = data.attributes.bandwidth(j);
        cfg.f0 = fc(j);
        clear tmp m;
        %% Beamforming
        [beams{j},ranges_out] = das_plane_wave_fd(cfg, sig_mf_dec);
        end
        timing_data(i) = toc;
        clear ping_data sig_mf_dec;
        if show_plots
            beam = [beams{1},beams{2},beams{3},beams{4},beams{5},beams{6}];
            sector_plot(data.attributes.scan_angle * pi/180,ranges_out, abs(beam), 20,"sector");
        end
     
    end
end
