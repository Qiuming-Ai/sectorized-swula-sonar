function timing_data = BaselineProcessFD(input_file, max_pings)
% Legacy desktop Dense workload. Timings are seconds; HDF5 loading is outside.
legacy_sonarsim_setup();
if nargin < 1 || isempty(input_file)
    paths = legacy_desktop_timing_inputs();
    input_file = paths.dense;
end
if nargin < 2, max_pings = Inf; end
[data, count] = legacy_load_timing_fixture(input_file, "Baseline", max_pings);
win = hamming(data.attributes.receive_array_num);
cfg.rx_xyz = data.attributes.receive_array_position;
cfg.c = data.attributes.sound_velocity;
cfg.t0 = data.attributes.sample_delay-data.attributes.pulse_duration/2;                      % 用你的实际 t_start
cfg.fs = data.attributes.sampling_frequency/data.attributes.decimate_factor;
cfg.angles = (-60:1:60) * pi/180;   % 扫描角（rad）
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
    %% Match Filter

    
    sig_mf = zeros(size(ping_data));
    for m = 1:data.attributes.receive_array_num
        tmp = conv(ping_data(:,m), data.attributes.match_filter_data, 'same');
        sig_mf(:,m) = tmp*win(m);
    end
    clear tmp ping_data m;
    
    %% TVG
    % To-do
    
    %% Beamforming
    [beam,ranges_out] = das_plane_wave_fd(cfg, sig_mf);
    timing_data(i) = toc;
    % sector_plot(cfg.angles,ranges_out, abs(beam), 20,"sector");
    % pause(0.1);
end
end
