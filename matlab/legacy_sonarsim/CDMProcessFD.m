function timing_data = CDMProcessFD(input_file, max_pings)
% Legacy desktop CDM workload. Timings are seconds; HDF5 loading is outside.
legacy_sonarsim_setup();
if nargin < 1 || isempty(input_file)
    paths = legacy_desktop_timing_inputs();
    input_file = paths.cdm;
end
if nargin < 2, max_pings = Inf; end
[data, count] = legacy_load_timing_fixture(input_file, "CDM", max_pings);
    win = hamming(data.attributes.receive_array_num);
    cfg.rx_xyz = data.attributes.receive_array_position;
    cfg.c = data.attributes.sound_velocity;
    cfg.t0 = data.attributes.sample_delay-data.attributes.pulse_duration/2;                      % 用你的实际 t_start
    cfg.fs = data.attributes.sampling_frequency/data.attributes.decimate_factor;
    % cfg.angles = data.attributes.scan_angle * pi/180;   % 扫描角（rad）
    cfg.BW = data.attributes.bandwidth;
    cfg.f0 = data.attributes.center_frequency;
    cfg.fd.nfft   =  2^nextpow2(length(data.pings{1}(:,1)));          % 或 64 / 256，越大分数时延越精细，但更慢
    cfg.fd.sign   = 1;
    cfg.cf.enable = false;
    cfg.cf.mode = 'cf';
    cfg.cf.gamma = 0.05;
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
        %% TVG
        % To-do
        index = data.attributes.scan_angle>=data.attributes.sector_div(j)&data.attributes.scan_angle<data.attributes.sector_div(j+1);
        cfg.angles =  data.attributes.scan_angle(index)/180*pi;
        clear tmp m;
        %% Beamforming
        [beams{j},ranges_out] = das_plane_wave_fd(cfg, sig_mf);
        end
        timing_data(i) = toc;
        clear ping_data;
        % beam = [beams{1},0.7*beams{2},1.5*beams{3},1.5*beams{4},0.7*beams{5},beams{6}];
        % sector_plot(data.attributes.scan_angle * pi/180,ranges_out, abs(beam), 20,"sector");
        % 
    end
end
