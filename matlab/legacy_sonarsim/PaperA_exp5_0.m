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


%% Baseline
probe.N                 = 128;             % Number of elements
pulse_len=5e-3; 
lightPos = [0, 0, 0];
rx_xyz = [linspace(-probe.N*probe.pitch/2,probe.N*probe.pitch/2,probe.N)',zeros(probe.N,2)];
% 绘制波束图（扫描角度0度）
[beamdb1,beam1, angles1] = plot_beam_pattern(rx_xyz, f0, ...
    'steer_angle', 0, ...
    'angle_range', [-85, 85], ...
    'angle_step', 0.1, ...
    'c', c0, ...
    'window', 'hamming', ...
    'plot_type', 'rect', ...
    'dyn_range', 120);
tw_db0 = beamdb1;
PSLR = calculate_PSLR(beamdb1, 0, angles1);
hold on;
%% CDM
Nrx = 32;             % Number of elements
Ntx = 12;

lightPos = [0, 0, 0];
rx_xyz = [linspace(-(Nrx-1)*lambda,(Nrx-1)*lambda,Nrx)',zeros(Nrx,2)];
tx_xyz = [linspace(-(Ntx-1)*lambda/4,(Ntx-1)*lambda/4,Ntx)',zeros(Ntx,2)];
angles_deg = -[-47.5,-26,-8.5,8.5,26,47.5]';              % 3 个方位角
fc   = 450e3;              % kHz 级别
BW   =100e3;

[exc_nt,exc_ori, delays_samp] = gen_multicode_cazac_excitation( ...
    tx_xyz, angles_deg, fc, BW, interp_factor*fs, c0,'amp',[1.2 1.2 1 1 1.2 1.2]);
 MF1 = conj(flipud(hilbert(exc_ori)));                  % FIR 形式
% pulse_len=length(MF(:,1))/(fs*interp_factor); 
% 假设你已有：
% X  : [Nt x M]   接收或发射端的多通道时间序列
% H  : [L  x Na]  每列一个匹配滤波器
% RX : [Nt x 3]   阵元坐标
% fs, c
[beamdb2,beam2, angles2] = plot_beam_pattern(rx_xyz, f0, ...
    'steer_angle', 0, ...
    'angle_range', [-85, 85], ...
    'angle_step', 0.1, ...
    'c', c0, ...
    'window', 'hamming', ...
    'plot_type', 'rect', ...
    'dyn_range', 60);
first_null = find_first_null(beamdb2,0,angles2);
index = angles2>0;
[BP_dB, theta_deg] = mf_beamscan_peakdir(exc_nt, MF1, tx_xyz, interp_factor*fs, c0);
 BP_dB = interp1(theta_deg,BP_dB',angles2);
%%
tx_bp = BP_dB(:,3).*(~index)+BP_dB(:,4).*(index);
theta_deg = angles2;
plot(theta_deg,tx_bp);
tw_db1 = tx_bp+beamdb2;

%% FDM
% 子带参数
fc   = [412e3; 458e3; 474e3; 490e3; 442e3; 426e3;];              % kHz 级别
BW   = [20e3; 20e3; 20e3;20e3;20e3;20e3];
pol  = ["up"; "up"; "up"; "up"; "up"; "up"];               % 上/下扫
pw   = repmat(pulse_len,[6 1]);           % 1.5 ms


[exc_nt,exc_ori, delays_samp] = gen_multisubband_lfm_excitation( ...
    tx_xyz, angles_deg, fc, BW, pol, pw, interp_factor*fs, c0);
 MF2 = conj(flipud(hilbert(exc_ori)));                  % FIR 形式


[BP_dB, theta_deg] = mf_beamscan_peakdir(exc_nt, MF2, tx_xyz, interp_factor*fs, c0);
tx_bp = interp1(theta_deg,BP_dB(3,:),angles2).*(~index)+interp1(theta_deg,BP_dB(4,:),angles2).*(index);
theta_deg = angles2;
plot(theta_deg,tx_bp);
tw_db2 = tx_bp+beamdb2;



pulse.fractional_bandwidth = 2/9;        % probe bandwidth [1]
pulse.center_frequency = f0;
BW=pulse.fractional_bandwidth*f0; 
impulse_response = 1;
t_p=0:1/fs:pulse_len;
excitation = chirp(t_p,f0-BW/2,pulse_len,f0+BW/2).*hamming(numel(t_p))';
one_way_ir = conv(impulse_response,excitation);
MF = conj(flipud(hilbert(excitation(:))));                  % FIR 形式
two_way_ir0 = conv(one_way_ir,MF);


two_way_ir1 = conv(conj(flipud(MF1(:,1))),MF1(:,1));
two_way_ir2 = conv(conj(flipud(MF2(:,1))),MF2(:,1));

two_way_ir1 = circshift(two_way_ir1,-480);

tw_db0_d = 20*log10(abs(two_way_ir0)/max(abs(two_way_ir0(:))));
tw_db1_d = resample(20*log10(abs(two_way_ir1)/max(abs(two_way_ir1(:)))),1,interp_factor);
tw_db2_d = resample(20*log10(abs(two_way_ir2)/max(abs(two_way_ir2(:)))),1,interp_factor);

l0 = length(tw_db0_d);l1 = length(tw_db1_d);l2 = length(tw_db2_d);
l = min([l0,l1,l2]);
tw_db0_d = tw_db0_d(1:l);
tw_db1_d = tw_db1_d(1:l);
tw_db2_d = tw_db2_d(1:l);
t_tw = [1:length(tw_db0_d)]/fs-pulse_len;
plot(t_tw,tw_db0_d);hold on;
plot(t_tw,tw_db1_d);
plot(t_tw,tw_db2_d);

