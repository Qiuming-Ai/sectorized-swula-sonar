%% 通带 LFM 信号脉冲压缩示例
% 本脚本演示线性调频（Linear Frequency Modulation, LFM）通带信号的
% 生成、目标回波仿真以及匹配滤波脉冲压缩处理流程。
% 场景采用典型水声参数：声速 1500 m/s、载频 450 kHz、带宽 20 kHz。
% 通过对比下变频与直接匹配滤波的结果，展示基带处理的重要性。

clear; clc; close all;

%% 参数设置
c = 1500;             % 水声波速 (m/s)
fc = 450e3;           % 载频 (Hz)
B = 20e3;             % 调频带宽 (Hz)
Tp = 10e-3;           % 脉冲宽度 (s)
Fs = 2e6;             % 采样率 (Hz)，需满足奈奎斯特采样

% 目标距离与回波延迟
R_target = 300;                     % 目标距离 (m)
tau = 2 * R_target / c;             % 往返传播时延 (s)

% 观测时间长度需覆盖发射脉冲与回波
guard_time = Tp;                    % 额外留出的观测余量
T_total = Tp + tau + guard_time;    % 总观测时间
N = ceil(T_total * Fs);             % 采样点数
t = (0:N-1) / Fs;                   % 时间轴
delay_samples = round(tau * Fs);    % 延迟对应的采样点

SNR_dB = 15;                        % 回波信噪比 (dB)

%% 生成通带 LFM 发射信号
k = B / Tp;                         % 调频斜率 (Hz/s)
t0 = t - Tp / 2;                    % 以脉冲中心对齐的时间轴
window = double(abs(t0) <= Tp / 2); % 矩形窗定义脉冲包络
lfm_baseband = window .* exp(1j * pi * k .* t0.^2); % 复包络

% 升余弦窗可替换矩形窗 (可选)
% beta = 0.25;
% window = double(abs(t0) <= Tp / 2) .* (0.5 * (1 + cos(pi * t0 / (Tp/2))));

tx_passband = real(lfm_baseband .* exp(1j * 2 * pi * fc .* t0)); % 通带实信号

%% 构造目标回波并添加噪声
rx_passband = zeros(1, length(tx_passband));
if delay_samples < length(tx_passband)
    rx_passband(delay_samples+1:end) = tx_passband(1:end-delay_samples);
end

% 添加白噪声
rx_passband = awgn(rx_passband, SNR_dB, 'measured');

%% 下变频至基带
analytic_rx = hilbert(rx_passband);                     % 构造解析信号
rx_baseband = analytic_rx .* exp(-1j * 2 * pi * fc .* t0);

%% 下采样降低基带采样率
target_Fs_ds = 40e3;                                    % 目标基带采样率，接近信号带宽
decim_factor = max(1, round(Fs / target_Fs_ds));
Fs_ds = Fs / decim_factor;

rx_baseband_ds = decimate(real(rx_baseband), decim_factor) + ...
                 1j * decimate(imag(rx_baseband), decim_factor);
lfm_baseband_ds = decimate(real(lfm_baseband), decim_factor) + ...
                  1j * decimate(imag(lfm_baseband), decim_factor);

%% 匹配滤波脉冲压缩（基带 + 下采样）
matched_filter = conj(fliplr(lfm_baseband_ds));         % 匹配滤波器
pulse_compressed = conv(rx_baseband_ds, matched_filter);% 脉压输出

% 不进行下变频的直接脉压
matched_filter_passband = fliplr(tx_passband);          % 实信号匹配滤波器
pulse_compressed_passband = conv(rx_passband, matched_filter_passband);

% 建立对应的时间轴
tc = (0:length(pulse_compressed)-1) / Fs_ds - tau;
tc_real = (0:length(pulse_compressed_passband)-1) / Fs - tau;
t_ds = (0:length(rx_baseband_ds)-1) / Fs_ds;

%% 绘图展示
figure('Name', 'LFM 脉冲压缩', 'Color', 'w', 'Position', [100 100 980 680]);

ax1 = subplot(4,1,1);
plot(t * 1e3, tx_passband);
xlabel('时间 / ms');
ylabel('幅度');
title('通带 LFM 发射波形');
grid on;

ax2 = subplot(4,1,2);
plot(t_ds * 1e3, abs(rx_baseband_ds));
xlabel('时间 / ms');
ylabel('幅度');
title(sprintf('下变频并下采样后的回波包络 (R = %d)', decim_factor));
grid on;

ax3 = subplot(4,1,3);
plot(tc * 1e3, abs(pulse_compressed) / max(abs(pulse_compressed)));
xlabel('时间 / ms');
ylabel('归一化幅度');
title('匹配滤波脉冲压缩结果');
grid on;

ax4 = subplot(4,1,4);
plot(tc_real * 1e3, abs(pulse_compressed_passband) / max(abs(pulse_compressed_passband)));
xlabel('时间 / ms');
ylabel('归一化幅度');
title('不下变频直接脉压输出（保留载频震荡）');
grid on;

% 调整子图宽度，为右侧注释预留空间
axes_array = [ax1, ax2, ax3, ax4];
for ax = axes_array
    pos = get(ax, 'Position');
    pos(1) = 0.08;
    pos(3) = 0.68;  % 缩短宽度
    set(ax, 'Position', pos);
end

% 在右侧标注采样率、信号频率及实/复信号
info_x = 0.78;
info_w = 0.2;

annotation('textbox', [info_x 0.74 info_w 0.12], ...
    'String', sprintf(['采样率: %.1f MHz\n' ...
                       '信号类型: 实信号\n' ...
                       '中心频率: %.1f kHz'], Fs/1e6, fc/1e3), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 10);

annotation('textbox', [info_x 0.53 info_w 0.12], ...
    'String', sprintf(['采样率: %.1f kHz\n' ...
                       '信号类型: 复信号\n' ...
                       '频率范围: ±%.1f kHz'], Fs_ds/1e3, B/2e3), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 10);

annotation('textbox', [info_x 0.32 info_w 0.12], ...
    'String', sprintf(['采样率: %.1f kHz\n' ...
                       '信号类型: 复信号\n' ...
                       '输出频率: 基带 (0 Hz)\n' ...
                       '不下变频: 保留载频振荡'], Fs_ds/1e3), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 10);

annotation('textbox', [info_x 0.11 info_w 0.12], ...
    'String', sprintf(['采样率: %.1f MHz\n' ...
                       '信号类型: 实信号\n' ...
                       '输出为实值振荡'], Fs/1e6), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 10);

sgtitle('通带 LFM 信号脉压仿真');

%% 计算并输出主瓣宽度（近似）
[~, peak_idx] = max(abs(pulse_compressed));
half_max = max(abs(pulse_compressed)) / sqrt(2);
indices = find(abs(pulse_compressed) >= half_max);
if numel(indices) >= 2
    mainlobe_width = (indices(end) - indices(1)) / Fs_ds;
    fprintf('脉压主瓣宽度约为 %.3f 微秒。\n', mainlobe_width * 1e6);
else
    fprintf('无法估计主瓣宽度。\n');
end

fprintf('目标距离：%.1f km，对应往返时延 %.2f 微秒。\n', R_target / 1e3, tau * 1e6);
fprintf('基带处理下采样因子 R = %d，等效采样率 %.1f kHz。\n', decim_factor, Fs_ds/1e3);

% 对比说明
fprintf(['下变频并下采样后的脉压输出保留基带包络，便于主瓣宽度和峰值读出；\n' ...
         '直接对通带信号匹配滤波会出现载频级联振荡，需要额外取包络或整形。\n']);


