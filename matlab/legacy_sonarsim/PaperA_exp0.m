% clear all; close all;
legacy_sonarsim_setup();
%% 閸╃儤婀扮敮鎼佸櫤
c0 = 1500;
fs = 2e6;
interp_factor = 10;
dt = 1/fs;
f0 = 0.45e+06;      % Transducer center frequency [Hz]

lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 15e-3;            % Height of element [m]
probe.pitch             = 0.5*lambda;        % probe.pitch [m]
kerf                    = 0.05*lambda;        % gap between elements [m]
probe.element_width     = probe.pitch-kerf;% Width of element [m]
% lens_el                 = 20e-3;           % position of the elevation focus
probe.N                 = 128;             % Number of elements
pulse_len=2e-3; 
lightPos = [0, 0, 0];
rx_xyz = [linspace(-probe.N*probe.pitch/2,probe.N*probe.pitch/2,probe.N)',zeros(probe.N,2)];
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

%% 濠碘偓閸旇精鍓﹂崘?
pulse.fractional_bandwidth = 2/9;        % probe bandwidth [1]
pulse.center_frequency = f0;
BW=pulse.fractional_bandwidth*f0; 
impulse_response = 1;
t_p=0:1/fs:pulse_len;
excitation = chirp(t_p,f0-BW/2,pulse_len,f0+BW/2).*hamming(numel(t_p))';
one_way_ir = conv(impulse_response,excitation);
MF = conj(flipud(excitation(:)));                  % FIR 瑜般垹绱?
two_way_ir = conv(one_way_ir,MF);
lag = length(two_way_ir)/2+1;   
%%
opts = struct('round','round','precision','single');
opts.array.use_resample = true;                   % 閼?true 鐏忔繆鐦拫鍐暏 MATLAB resample

% opts.array.density_scale_tx = 1;   % 闂傜绐涜劤2閿涘牆鐦戞惔锔惧閸戝繐宕愰敍?
% opts.array.density_scale_rx = 1;   % 闂傜绐涜劤3
% opts.array.max_sub_per_elem_tx = 128; % 濮ｅ繋閲淭X闂冮潧鍘撻惃鍕珓閹风喎鐡欓梼鍏哥瑐闂?
% opts.array.max_sub_per_elem_rx = 256; % 濮ｅ繋閲淩X闂冮潧鍘撴稉濠囨
% opts.array.rx.shape  = 'disk';
% opts.array.rx.size   = [probe.width/2 probe.height];
% opts.array.rx.steer = [45*ones(probe.N,1) zeros(probe.N,2)];
[v, t] = sim_rx_from_scatterers_simple(point_position,point_amplitudes,lightPos,rx_xyz,f0,c0,fs,excitation,opts);

%% Match Filter
sig_mf = zeros(size(v));
win = hamming(probe.N);
for m = 1:probe.N
    tmp = conv(v(:,m), MF, 'same');
    sig_mf(:,m) = tmp*win(m);
end
% TVG
tvg = tvg_factory(t,'c',c0,'mode','custom','freq_kHz',f0/1e3,'G0_dB',0,'K',40);
% 閸嬪洩顔曢柌鍥ㄧ壉閻?fs=500 kHz閿涘矁顔囪ぐ?100 ms閿?
tau = (0:length(sig_mf)-1)/fs;  % 锜挎潪?
G_amp = tvg.Gamp(tau)';  % 閻劋绨獮鍛/閸栧懐绮剁悰銉ヤ缉
for m = 1:probe.N
    sig_mf(:,m) = sig_mf(:,m).*G_amp;
end

tic
% %% Coherent Plane Wave DAS
% 濮濄儵顎?: 閹笛嗩攽闁插洦鐗遍悳鍥祮閹?
if interp_factor>1
converter = dsp.SampleRateConverter(...
    'InputSampleRate', fs, ...
    'OutputSampleRate', interp_factor*fs, ...
    'Bandwidth', 0.45 * fs, ... % 闁艾鐢敮锕€顔?
    'StopbandAttenuation', 80);          % 闂冭鐢悰鏉垮櫤 (dB)

sig_mf = converter(sig_mf);
end
%%
data = hilbert(sig_mf./max(sig_mf(:)));
rx = data;
cfg = struct;
cfg.interp_factor = interp_factor;
cfg.pulse_len = pulse_len;
cfg.impulse_response = impulse_response;
cfg.excitation = excitation;
cfg.rx_xyz = rx_xyz;
cfg.c = c0;
cfg.t0 = t-pulse_len/1.992;                      % 閻劋缍橀惃鍕杽闂?t_start
cfg.fs = interp_factor*fs;
cfg.is_baseband = false;    % 婵″倹鐏?X 閺勵垰鐔€鐢毌Q
cfg.fc = f0;            % 900kHz 缁€杞扮伐
cfg.nfft = 2^nextpow2(size(rx,1));    % 閹?1,048,576 缁?
cfg.winLen = 4096;
cfg.hopLen = 1024;
cfg.numSnap = 96;          % 绾喕绻?#韫囶偅濯?>= M
cfg.diagLoad = 1e-2;
cfg.tx_xyz = lightPos; % 閸欐垵鐨犻崷銊╂█閸掓ぞ鑵?
% cfg.rx_win = hann(128);                % 鐎涙柨绶為崝鐘虫綀閿涘矂妾烽弮浣烘憵
cfg.bw = BW;
cfg.angles = (-60:0.1:60) * pi/180;   % 閹殿偅寮跨憴鎺炵礄rad閿?
cfg.ranges = (4.5:0.001:6).';         % 閹殿偅寮垮ǎ鍗炲閿涘潰閿?
% 閸嬪洩顔?rx 娑?Nsamp鑴?4 閻ㄥ嫬甯慨瀣槻閸婂吋鏆熼幑?
cfg.cf.enable = false;
cfg.cf.mode   = 'cf';     % 閹?'pcf'
cfg.cf.gamma  = 1.2;      % >1 閺囨潙宸遍梻銊﹀付
cfg.cf.eps    = 1e-12;

cfg.sva.enable  = false;                         % 瀵偓/閸?SVA
cfg.sva.windows = {'rect','hann','hamming', 'kaiser:6'};  % 閼奉亜鐣炬稊澶岀崶閺?
cfg.sva.ref_idx = 1;                            % 閻╅晲缍呴崣鍌濃偓鍐崶閿涘牏鐓╄ぐ顫礆
cfg.sva.norm    = 'mean';                       % 濮ｅ繒鐛ラ崸鍥р偓鐓庣秺娑撯偓閿涘牊甯归懡鎰剁礆
img = das_plane_wave(cfg, rx);
toc
% 閸欘垵顫嬮崠鏍电礄瀵搫瀹砫B閿?
% 閸嬪洩顔曟担鐘插嚒閻劌澧犻棃銏㈡畱 DAS 瀵版鍩?img (Nr鑴砃a)
sector_plot(cfg.angles, cfg.ranges, img, 30,'rect');
title('Plane-wave DAS (integer delays)');hold on;
scatter3(point_position(:,1),point_position(:,3),point_position(:,2),20,'r','filled');
% 在指定距离处绘制“角度-能量”图，并用红色竖线标识目标角度
target_range_m = 5.12; % 修改这里即可切换目标距离（单位：m）
target_ranges_m = sqrt(point_position(:,1).^2 + point_position(:,3).^2);
target_angles_deg = atan2d(point_position(:,1), point_position(:,3));
range_gate_m = max(0.003, median(diff(cfg.ranges))); % 默认约 ±3 mm 门限
idx_target = abs(target_ranges_m - target_range_m) <= range_gate_m;
plot_angle_energy_at_range(cfg.angles, cfg.ranges, img, target_range_m, 40, target_angles_deg(idx_target));

function plot_angle_energy_at_range(angles_rad, ranges_m, img, target_range_m, dyn_range_db, target_angles_deg)
% 在指定距离处提取不同角度的能量并可视化（dB）
if nargin < 5 || isempty(dyn_range_db)
    dyn_range_db = 40;
end
if nargin < 6
    target_angles_deg = [];
end

if isempty(img) || isempty(angles_rad) || isempty(ranges_m)
    error("angles/ranges/img 不能为空。");
end

[~, ridx] = min(abs(ranges_m - target_range_m));
actual_range_m = ranges_m(ridx);

% 提取该距离处随角度变化的能量（强度）
energy = abs(img(ridx, :)).^2;
energy = energy ./ max(energy + eps);
energy_db = 10 * log10(energy + eps);
energy_db = max(energy_db, -abs(dyn_range_db));

figure;
plot(angles_rad * 180 / pi, energy_db, 'LineWidth', 1.5);
grid on;
xlabel('Angle (deg)');
ylabel('Energy (dB)');
title(sprintf('Energy vs Angle at range = %.3f m (nearest bin %.3f m)', target_range_m, actual_range_m));
xlim([min(angles_rad) max(angles_rad)] * 180 / pi);
ylim([-abs(dyn_range_db) 0]);
hold on;

% 用红色竖线标注目标角度（若有）
for k = 1:numel(target_angles_deg)
    xline(target_angles_deg(k), 'r--', 'LineWidth', 1.2);
end
end
