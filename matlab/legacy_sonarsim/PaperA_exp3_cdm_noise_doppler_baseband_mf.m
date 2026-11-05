% clear all; close all;
legacy_sonarsim_setup();
% load FIR_Filter.mat
%% 基本常量
c0 = 1500;
fs = 2e6;
interp_factor = 8;
dt = 1/fs;
f0 = 0.45e+06;      % Transducer center frequency [Hz]
BW=100e3;
fmax = f0+BW/2;
fmin = f0-BW/2;
lambda                  = c0/fmax;           % Wavelength [m]

Nrx = 32;             % Number of elements
Ntx = 12;

lightPos = [0, 0, 0];
rx_xyz = [linspace(-(Nrx-1)*lambda,(Nrx-1)*lambda,Nrx)',zeros(Nrx,2)];
tx_xyz = [linspace(-(Ntx-1)*lambda/4,(Ntx-1)*lambda/4,Ntx)',zeros(Ntx,2)];
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
% 子带参数
angles_deg = -[-47.5,-26,-8.5,8.5,26,47.5]';              % 3 个方位角
fc   = 450e3;              % kHz 级别
BW   =100e3;

[exc_nt,exc_ori, delays_samp] = gen_multicode_cazac_excitation( ...
    tx_xyz, angles_deg, fc, BW, interp_factor*fs, c0,'amp',[1.2 1.2 1 1 1.2 1.2]);
 MF = conj(flipud(hilbert(exc_ori)));                  % FIR 形式
pulse_len=length(MF(:,1))/(fs*interp_factor); 
% 假设你已有：
% X  : [Nt x M]   接收或发射端的多通道时间序列
% H  : [L  x Na]  每列一个匹配滤波器
% RX : [Nt x 3]   阵元坐标
% fs, c

[BP_dB, theta_deg] = mf_beamscan_peakdir(exc_nt, MF, tx_xyz, interp_factor*fs, c0);
% % 
% % 画出第 a 个子带的指向性：
figure; 
for i = 1:6
a = i;
plot(theta_deg, BP_dB(a,:), 'LineWidth', 1.5);hold on;
end
xlabel('Angle (deg)'); ylabel('Normalized magnitude (dB)');
title(sprintf('Subband %d beampattern', a)); grid on;legend({'Sub-Beam 1','Sub-Beam 2','Sub-Beam 3','Sub-Beam 4','Sub-Beam 5','Sub-Beam 6'})
ylim([-60, 3]);


%%
opts = struct('round','round','precision','single');
opts.array.use_resample = true;                   % 若 true 尝试调用 MATLAB resample

[v, t] = sim_rx_from_scatterers_perTX(point_position,point_amplitudes,tx_xyz,rx_xyz,f0,c0,interp_factor*fs,exc_nt,opts);
v = (v-min(v(:)))./(max(v(:))-min(v(:)));
v = (v-0.5)*2;
MF = resample(MF,1,interp_factor);
v =  resample(v,1,interp_factor);
figure;
win = hamming(Nrx);
filt = SignalFilter('SampleRate',fs, ...
    'Passband', BW/2, ...
    'StopbandAtten', 60, ...
    'Type', 'lowpass', ...
    'SidebandMode', 'single', ...
    'OutputSameLength', true);

MF_mix = MF.*(exp(1j*2*pi*fc*[0:length(MF(:,1))-1]/(fs))');
angles_div{1} = [-60:0.1:-34-0.1];
angles_div{2} = [-34:0.1:-17-0.1];
angles_div{3} = [-17:0.1:0-0.1];
angles_div{4} = [0:0.1:17-0.1];
angles_div{5} = [17:0.1:34-0.1];
angles_div{6} = [34:0.1:60];

%%
for k = 1:length(snr_level)
    for j = 1:length(velocity)
        % v_calc = resample(v,1,interp_factor);
        v_calc = resample(v,c0+2*velocity(j),c0);
        v_calc_1 = awgn(v_calc,snr_level(k),'measured');
        v_calc_2 = awgn(v_calc,snr_level(k),'measured');
        v_calc_3 = awgn(v_calc,snr_level(k),'measured');
        v_calc_4 = awgn(v_calc,snr_level(k),'measured');
        v_calc = v_calc_1+v_calc_2+v_calc_3+v_calc_4;
        beam = [];
        % tic
        %% 下变频
        v_mix = v_calc.*repmat(exp(1j*2*pi*fc*[0:length(v_calc(:,1))-1]/fs)',[1 Nrx]);
        for i = 1:6

            %% Match Filter
            sig_mf = zeros(size(v_mix));
            for m = 1:Nrx
                tmp = conv(v_mix(:,m), MF_mix(:,i), 'same');
                sig_mf(:,m) = tmp*win(m);
            end
                        %% FIR滤波器
            y_fir = filt.process(sig_mf);
            % MF_fir = filt.process(MF_mix(:,i));
            %% Decimation
            y_deci = resample(y_fir,1,16);
            % MF_deci =  resample(MF_fir,1,16);
            fs_deci = fs/16;


            %%
            data = y_deci./max(abs(y_deci(:)));
            rx = data;
            cfg = struct;
            cfg.rx_xyz = rx_xyz;
            cfg.c = c0;
            cfg.t0 = t-pulse_len/2;                      % 用你的实际 t_start
            cfg.fs = fs_deci;
            cfg.angles = angles_div{i} * pi/180;   % 扫描角（rad）
            cfg.ranges = (4.5:0.01:6).';         % 扫描深度（m）
            cfg.BW = 120e3;
            cfg.f0 = fc;
            % cfg.w = hamming(probe.N);
            % 假设 rx 为 Nsamp×64 的原始复值数据
            cfg.fd.nfft   = 1024;          % 或 64 / 256，越大分数时延越精细，但更慢
            cfg.fd.sign   = 1;
            cfg.cf.enable = true;
            cfg.cf.mode = 'cf';
            cfg.cf.gamma = 0.4;
            [beam{i},ranges_out] = das_plane_wave_fd(cfg, rx);
        end
        beams = [beam{1},0.7*beam{2},1.6*beam{3},1.6*beam{4},0.7*beam{5},beam{6}];
        % toc
        cfg.angles = [-60:0.1:60]/180*pi;
        % 可视化（强度dB）
        % 假设你已用前面的 DAS 得到 img (Nr×Na)
        ax = subplot(length(snr_level),length(velocity),(j-1)*3+k);
        % set(ax,'Parent',h(k))
        sector_plot(cfg.angles, ranges_out, beams, 20,"sector");
        % title('Plane-wave DAS (integer delays)');
        hold on;
        % scatter3(asin(point_position(:,1)./sqrt(point_position(:,1).^2+point_position(:,3).^2))/pi*180,sqrt(point_position(:,1).^2+point_position(:,3).^2),point_position(:,2),8,'r','filled');
        % xlim([-2 2]),ylim([4.5 6])
        ylim([4.8 5.7])
    end
end



