clear all; close all;
legacy_sonarsim_setup();
%% 基本常量
c0 = 1500;
fs = 2e6;
interp_factor = 8;
dt = 1/fs;
f0 = 0.45e+06;      % Transducer center frequency [Hz]

lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 15e-3;            % Height of element [m]
probe.pitch             = 0.5*lambda;        % probe.pitch [m]
kerf                    = 0.05*lambda;        % gap between elements [m]
probe.element_width     = probe.pitch-kerf;% Width of element [m]
% lens_el                 = 20e-3;           % position of the elevation focus
probe.N                 = 128;             % Number of elements
pulse_len=5e-3; 
lightPos = [0, 0, 0];
rx_xyz = [linspace(-probe.N*probe.pitch/2,probe.N*probe.pitch/2,probe.N)',zeros(probe.N,2)];
velocity = [0];
snr_level = [10];
%% Define phantom
% Define some points in a phantom for the simulation
point_position(1,:) = [0 0 498e-2];
% point_position(end+1,:) = [0 0 500e-2];
% point_position(end+1,:) = [0 0 504e-2];
% point_position(end+1,:) = [0 0 512e-2];
% point_position(end+1,:) = [0 0 528e-2];
% point_position(end+1,:) = [0 0 550e-2];
% point_position(end+1,:) = [512e-2*sind(2) 0 512e-2*cosd(2)];
% point_position(end+1,:) = [512e-2*sind(-2) 0 512e-2*cosd(-2)];
% point_position(end+1,:) = [512e-2*sind(6) 0 512e-2*cosd(6)];
% point_position(end+1,:) = [512e-2*sind(-6) 0 512e-2*cosd(-6)];
% point_position(end+1,:) = [512e-2*sind(12) 0 512e-2*cosd(12)];
% point_position(end+1,:) = [512e-2*sind(-12) 0 512e-2*cosd(-12)];
% point_position(end+1,:) = [512e-2*sind(20) 0 512e-2*cosd(20)];
% point_position(end+1,:) = [512e-2*sind(-20) 0 512e-2*cosd(-20)];
% point_position(end+1,:) = [512e-2*sind(36) 0 512e-2*cosd(36)];
% point_position(end+1,:) = [512e-2*sind(-36) 0 512e-2*cosd(-36)];
% point_position(end+1,:) = [512e-2*sind(42) 0 512e-2*cosd(42)];
% point_position(end+1,:) = [512e-2*sind(-42) 0 512e-2*cosd(-42)];
% point_position(end+1,:) = [512e-2*sind(55) 0 512e-2*cosd(55)];
% point_position(end+1,:) = [512e-2*sind(-55) 0 512e-2*cosd(-55)];
% Set point amplitudes
point_amplitudes = ones(size(point_position,1),1);

%% 激励脉冲
pulse.fractional_bandwidth = 2/9;        % probe bandwidth [1]
pulse.center_frequency = f0;
BW=pulse.fractional_bandwidth*f0; 
impulse_response = 1;
t_p=0:1/fs:pulse_len;
excitation = chirp(t_p,f0-BW/2,pulse_len,f0+BW/2).*hamming(numel(t_p))';
one_way_ir = conv(impulse_response,excitation);
MF = conj(flipud(hilbert(excitation(:))));                  % FIR 形式
two_way_ir = conv(one_way_ir,MF);
lag = length(two_way_ir)/2+1;   
%%
opts = struct('round','round','precision','single');
opts.array.use_resample = true;                   % 若 true 尝试调用 MATLAB resample

% opts.array.density_scale_tx = 1;   % 间距×2（密度约减半）
% opts.array.density_scale_rx = 1;   % 间距×3
% opts.array.max_sub_per_elem_tx = 128; % 每个TX阵元的虚拟子阵上限
% opts.array.max_sub_per_elem_rx = 256; % 每个RX阵元上限
% opts.array.rx.shape  = 'disk';
% opts.array.rx.size   = [probe.width/2 probe.height];
% opts.array.rx.steer = [45*ones(probe.N,1) zeros(probe.N,2)];
[v, t] = sim_rx_from_scatterers_simple(point_position,point_amplitudes,lightPos,rx_xyz,f0,c0,fs,excitation,opts);
v = (v-min(v(:)))./(max(v(:))-min(v(:)));
v = (v-0.5)*2;
figure;
win = hamming(probe.N);
filt = SignalFilter('SampleRate',fs, ...
    'Passband', BW/2, ...
    'StopbandAtten', 60, ...
    'Type', 'lowpass', ...
    'SidebandMode', 'single', ...
    'OutputSameLength', true);
MF_mix = MF.*(exp(1j*2*pi*f0*[0:length(MF)-1]/fs)');
for k = 1:length(snr_level)
    for j = 1:length(velocity)
        v_calc = resample(v,c0+2*velocity(j),c0);
        v_calc = awgn(v_calc,snr_level(k),'measured');
        %% 下变频
        v_mix = v_calc.*repmat(exp(1j*2*pi*f0*[0:length(v_calc(:,1))-1]/fs)',[1 probe.N]);
        %% FIR滤波器
        y_fir = filt.process(v_mix);
        MF_fir = filt.process(MF_mix);
        %% Decimation
        y_deci = resample(y_fir,1,2);
        MF_deci =  resample(MF_fir,1,2);
        % for i = 1:probe.N
        % y_deci(:,i) = decimate(y_fir(:,i),8);
        % end
        % MF_deci = decimate(MF_fir,8);
        fs_deci = fs/2;
        %% Match Filter
        sig_mf = zeros(size(y_deci));
        tic
        for m = 1:probe.N
            tmp = conv(y_deci(:,m), MF_deci, 'same');
            sig_mf(:,m) = tmp*win(m);
        end
        % toc
        % TVG
        tvg = tvg_factory(t,'c',c0,'mode','custom','freq_kHz',f0/1e3,'G0_dB',0,'K',40);
        % 假设采样率 fs=500 kHz，记录 100 ms：
        tau = (0:length(sig_mf)-1)/fs_deci;  % τ轴
        G_amp = tvg.Gamp(tau)';  % 用于幅度/包络补偿
        for m = 1:probe.N
            sig_mf(:,m) = sig_mf(:,m).*G_amp;
        end

        data = sig_mf./max(abs(sig_mf(:)));
        rx = data;
        cfg = struct;
        cfg.rx_xyz = rx_xyz;
        cfg.c = c0;
        cfg.t0 = t-pulse_len/2;                      % 用你的实际 t_start
        cfg.fs = fs_deci;
        cfg.angles = (-60:0.1:60) * pi/180;   % 扫描角（rad）
        cfg.ranges = (4.5:0.01:6).';         % 扫描深度（m）
        cfg.BW = 200e3;
        cfg.f0 = 450e3;
        % cfg.w = hamming(probe.N);
        % 假设 rx 为 Nsamp×64 的原始复值数据
        cfg.fd.nfft   = 1024;          % 或 64 / 256，越大分数时延越精细，但更慢
        cfg.fd.sign   = 1;
        cfg.cf.enable = true;
        cfg.cf.mode = 'cf';
        cfg.cf.gamma = 0.1;
        [beam,ranges_out] = das_plane_wave_fd(cfg, rx);
                % toc
        % 可视化（强度dB）
        % 假设你已用前面的 DAS 得到 img (Nr×Na)
        ax = subplot(length(snr_level),length(velocity),(j-1)*3+k);
        % sector_plot(cfg.angles, cfg.ranges, img, 30,"rect");
        beam = compensate_range_shift(beam, 1/((cfg.fd.nfft/cfg.fs)/c0), c0, 0.035,'true');
        sector_plot(cfg.angles,ranges_out, abs(beam), 40,"sector");
        % title('Plane-wave DAS (integer delays)');
        hold on;
        % scatter3(point_position(:,1),point_position(:,3),point_position(:,2),8,'r','filled');
        % scatter3(asin(point_position(:,1)./sqrt(point_position(:,1).^2+point_position(:,3).^2))/pi*180,sqrt(point_position(:,1).^2+point_position(:,3).^2),point_position(:,2),12,'r','filled');
         % ylim([4.8 5.7]);
                 view(3)
        xlim([-1 1])
        ylim([4.9 5.05])
    end
end




