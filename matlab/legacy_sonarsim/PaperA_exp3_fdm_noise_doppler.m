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

lambda                  = c0/fmax;           % Wavelength [m]

Nrx                 = 32;             % Number of elements
Ntx = 12;
pulse_len=0.5e-3; 
lightPos = [0, 0, 0];
rx_xyz = [linspace(-(Nrx-1)*lambda,(Nrx-1)*lambda,Nrx)',zeros(Nrx,2)];
tx_xyz = [linspace(-(Ntx-1)*lambda/4,(Ntx-1)*lambda/4,Ntx)',zeros(Ntx,2)];
velocity = [1 5 10];
snr_level = [-5 -10 -15];
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
% 子带参数
angles_deg = -[-47.5,-26,-8.5,8.5,26,47.5]';              % 3 个方位角
fc   = [410e3; 458e3; 474e3; 490e3; 442e3; 426e3;];              % kHz 级别
BW   = [20e3; 20e3; 20e3;20e3;20e3;20e3];
pol  = ["up"; "up"; "up"; "up"; "up"; "up"];               % 上/下扫
pw   = repmat(pulse_len,[6 1]);           % 1.5 ms


[exc_nt,exc_ori, delays_samp] = gen_multisubband_lfm_excitation( ...
    tx_xyz, angles_deg, fc, BW, pol, pw, interp_factor*fs, c0);
 MF = conj(flipud(exc_ori));                  % FIR 形式

% 假设你已有：
% X  : [Nt x M]   接收或发射端的多通道时间序列
% H  : [L  x Na]  每列一个匹配滤波器
% RX : [Nt x 3]   阵元坐标
% fs, c

% [BP_dB, theta_deg] = mf_beamscan_peakdir(exc_nt, MF, tx_xyz, interp_factor*fs, c0);
% % 
% % 画出第 a 个子带的指向性：
% figure; 
% for i = 1:6
% a = i;
% plot(theta_deg, BP_dB(a,:), 'LineWidth', 1.5);hold on;
% end
% xlabel('Angle (deg)'); ylabel('Normalized magnitude (dB)');
% title(sprintf('Subband %d beampattern', a)); grid on;legend({'Sub-Beam 1','Sub-Beam 2','Sub-Beam 3','Sub-Beam 4','Sub-Beam 5','Sub-Beam 6'})
% ylim([-60, 3]);


%%
opts = struct('round','round','precision','single');
opts.array.use_resample = true;                   % 若 true 尝试调用 MATLAB resample

[v, t] = sim_rx_from_scatterers_perTX(point_position,point_amplitudes,tx_xyz,rx_xyz,f0,c0,interp_factor*fs,exc_nt,opts);
v = (v-min(v(:)))./(max(v(:))-min(v(:)));
v = (v-0.5)*2;
figure;
win = hamming(Nrx);
converter = dsp.SampleRateConverter(...
    'InputSampleRate', fs, ...
    'OutputSampleRate', interp_factor*fs, ...
    'Bandwidth', 0.45 * fs, ... % 通带带宽
    'StopbandAttenuation', 80);          % 阻带衰减 (dB)
v = hilbert(v);
MF = hilbert(MF);
for k = 1:3
    for j = 1:3
        v_calc = resample(v,1,interp_factor);
        v_calc = resample(v_calc,c0+2*velocity(j),c0);
        v_calc_1 = awgn(v_calc,snr_level(k),'measured');
        v_calc_2 = awgn(v_calc,snr_level(k),'measured');
        v_calc_3 = awgn(v_calc,snr_level(k),'measured');
        v_calc_4 = awgn(v_calc,snr_level(k),'measured');
        v_calc = v_calc_1+v_calc_2+v_calc_3+v_calc_4;
        MF_calc = resample(MF,1,interp_factor);
        angles_div{1} = [-60:0.1:-34-0.1];
        angles_div{2} = [-34:0.1:-17-0.1];
        angles_div{3} = [-17:0.1:0-0.1];
        angles_div{4} = [0:0.1:17-0.1];
        angles_div{5} = [17:0.1:34-0.1];
        angles_div{6} = [34:0.1:60];
        img = [];
        % tic
        for i = 1:6
            %% Match Filter
            sig_mf = zeros(size(v_calc));
            % for ii = 1:numel(filters)
            %     if filters(ii).id==i
            %         firData = filters(ii).b;
            %     end
            % end
            % tic
            for m = 1:Nrx
                % tmp = conv(v_calc(:,m), firData, 'same');
                % figure,
                % freqz(firData)
                % figure,
                % freqz(MF_calc(:,i))
                tmp = conv(v_calc(:,m), MF_calc(:,i), 'same');
                sig_mf(:,m) = tmp*win(m);
            end
            % echo("Match ")
            % toc
            % echo("\n")
            % %% Coherent Plane Wave DAS
            % 步骤4: 执行采样率转换
            if interp_factor>1
                sig_mf = converter(sig_mf);
            end
            %%
            % data = hilbert(sig_mf);

            rx = sig_mf;
            cfg = struct;
            cfg.interp_factor = interp_factor;
            cfg.pulse_len = pulse_len;
            cfg.rx_xyz = rx_xyz;
            cfg.c = c0;
            cfg.t0 = t-pulse_len/1.992;                      % 用你的实际 t_start
            cfg.fs = interp_factor*fs;
            cfg.fc = f0;            % 900kHz 示例
            cfg.tx_xyz = mean(tx_xyz); % 发射在阵列中
            % cfg.rx_win = hann(128);                % 孔径加权，降旁瓣
            cfg.bw = BW;
            cfg.angles = angles_div{i} * pi/180;   % 扫描角（rad）
            cfg.ranges = (4.8:0.001:5.7).';         % 扫描深度（m）
            tic
            img{i} = das_plane_wave(cfg, rx);
            % echo("DAS ")
            toc
            % echo("\n")
        end
        imgs = [img{1},img{2},img{3},img{4},img{5},img{6}];
        % toc
        cfg.angles = [-60:0.1:60]/180*pi;
        % 可视化（强度dB）
        % 假设你已用前面的 DAS 得到 img (Nr×Na)
        ax = subplot(3,3,(j-1)*3+k);
        % set(ax,'Parent',h(k))
        sector_plot(cfg.angles, cfg.ranges, imgs, 30,"rect");
        % title('Plane-wave DAS (integer delays)');
        hold on;
        scatter3(asin(point_position(:,1)./sqrt(point_position(:,1).^2+point_position(:,3).^2))/pi*180,sqrt(point_position(:,1).^2+point_position(:,3).^2),point_position(:,2),8,'r','filled');
        % xlim([-2 2]),ylim([4.5 6])
    end
end



