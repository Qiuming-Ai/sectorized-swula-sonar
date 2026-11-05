% clear all; close all;
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
snr_level = [-10 -20 -30];
%% Define phantom
% Define some points in a phantom for the simulation
point_position(1,:) = [0 0 498e-2];
point_position(end+1,:) = [0 0 500e-2];
point_position(end+1,:) = [0 0 504e-2];
point_position(end+1,:) = [0 0 512e-2];
point_position(end+1,:) = [0 0 528e-2];
point_position(end+1,:) = [0 0 550e-2];
point_position(end+1,:) = [512e-2*sind(2) 0 512e-2*cosd(2)];
point_position(end+1,:) = [512e-2*sind(-2) 0 512e-2*cosd(-2)];
point_position(end+1,:) = [512e-2*sind(6) 0 512e-2*cosd(6)];
point_position(end+1,:) = [512e-2*sind(-6) 0 512e-2*cosd(-6)];
point_position(end+1,:) = [512e-2*sind(12) 0 512e-2*cosd(12)];
point_position(end+1,:) = [512e-2*sind(-12) 0 512e-2*cosd(-12)];
point_position(end+1,:) = [512e-2*sind(20) 0 512e-2*cosd(20)];
point_position(end+1,:) = [512e-2*sind(-20) 0 512e-2*cosd(-20)];
point_position(end+1,:) = [512e-2*sind(36) 0 512e-2*cosd(36)];
point_position(end+1,:) = [512e-2*sind(-36) 0 512e-2*cosd(-36)];
point_position(end+1,:) = [512e-2*sind(42) 0 512e-2*cosd(42)];
point_position(end+1,:) = [512e-2*sind(-42) 0 512e-2*cosd(-42)];
point_position(end+1,:) = [512e-2*sind(55) 0 512e-2*cosd(55)];
point_position(end+1,:) = [512e-2*sind(-55) 0 512e-2*cosd(-55)];
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
MF = conj(flipud(excitation(:)));                  % FIR 形式
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
converter = dsp.SampleRateConverter(...
    'InputSampleRate', fs, ...
    'OutputSampleRate', interp_factor*fs, ...
    'Bandwidth', 0.45 * fs, ... % 通带带宽
    'StopbandAttenuation', 80);          % 阻带衰减 (dB)
v = hilbert(v);
MF = hilbert(MF);
for k = 1:length(snr_level)
    for j = 1:length(velocity)
        v_calc = resample(v,c0+2*velocity(j),c0);
        v_calc = awgn(v_calc,snr_level(k),'measured');
        %% Match Filter
        sig_mf = zeros(size(v_calc));
        tic
        for m = 1:probe.N
            tmp = conv(v_calc(:,m), MF, 'same');
            sig_mf(:,m) = tmp*win(m);
        end
        % toc
        % TVG
        tvg = tvg_factory(t,'c',c0,'mode','custom','freq_kHz',f0/1e3,'G0_dB',0,'K',40);
        % 假设采样率 fs=500 kHz，记录 100 ms：
        tau = (0:length(sig_mf)-1)/fs;  % τ轴
        G_amp = tvg.Gamp(tau)';  % 用于幅度/包络补偿
        for m = 1:probe.N
            sig_mf(:,m) = sig_mf(:,m).*G_amp;
        end

        % %% Coherent Plane Wave DAS
        % 步骤4: 执行采样率转换
        % tic
        if interp_factor>1
            sig_mf = converter(sig_mf);
        end
        %%
        % data = hilbert(sig_mf./max(sig_mf(:)));
        rx = sig_mf;
        cfg = struct;
        cfg.interp_factor = interp_factor;
        cfg.pulse_len = pulse_len;
        cfg.impulse_response = impulse_response;
        cfg.excitation = excitation;
        cfg.rx_xyz = rx_xyz;
        cfg.c = c0;
        cfg.t0 = t-pulse_len/1.992;                      % 用你的实际 t_start
        cfg.fs = interp_factor*fs;
        cfg.is_baseband = false;    % 如果 X 是基带IQ
        cfg.fc = f0;            % 900kHz 示例
        cfg.nfft = 2^nextpow2(size(rx,1));    % 或 1,048,576 等
        cfg.winLen = 4096;
        cfg.hopLen = 1024;
        cfg.numSnap = 96;          % 确保 #快拍 >= M
        cfg.diagLoad = 1e-2;
        cfg.tx_xyz = lightPos; % 发射在阵列中
        % cfg.rx_win = hann(128);                % 孔径加权，降旁瓣
        cfg.bw = BW;
        cfg.angles = (-60:0.1:60) * pi/180;   % 扫描角（rad）
        cfg.ranges = (4.8:0.001:5.7).';         % 扫描深度（m）
        % 假设 rx 为 Nsamp×64 的原始复值数据
        cfg.cf.enable = true;
        cfg.cf.mode   = 'cf';     % 或 'pcf'
        cfg.cf.gamma  = 0.2;      % >1 更强门控
        cfg.cf.eps    = 1e-12;

        cfg.sva.enable  = false;                         % 开/关 SVA
        cfg.sva.windows = {'rect','hann','hamming', 'kaiser:6'};  % 自定义窗族
        cfg.sva.ref_idx = 1;                            % 相位参考窗（矩形）
        cfg.sva.norm    = 'mean';                    % 每窗均值归一（推荐）
        tic  
        img = das_plane_wave(cfg, rx);
        toc
        % toc
        % 可视化（强度dB）
        % 假设你已用前面的 DAS 得到 img (Nr×Na)
        ax = subplot(length(snr_level),length(velocity),(j-1)*3+k);
        sector_plot(cfg.angles, cfg.ranges, img, 20,"rect");
        % title('Plane-wave DAS (integer delays)');
        hold on;
        % scatter3(point_position(:,1),point_position(:,3),point_position(:,2),8,'r','filled');
        scatter3(asin(point_position(:,1)./sqrt(point_position(:,1).^2+point_position(:,3).^2))/pi*180,sqrt(point_position(:,1).^2+point_position(:,3).^2),point_position(:,2),8,'r','filled');
    end
end




