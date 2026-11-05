function stats = fd_fir_delay_selftest(mu, L, method, beta, Fs, B, N, doplot)
% fd_fir_delay_selftest
% 单文件自测：用 FIR 分数时延实现 y[n]≈x[n-mu]，估计实际时延并评估误差。
%
% USAGE:
%   stats = fd_fir_delay_selftest();                  % 跑默认示例并出图
%   stats = fd_fir_delay_selftest(0.35,65,'sinc',6,1e6,2e5,65536,true);
%
% INPUTS:
%   mu      : 目标总时延（样点，可非整数、可正负）
%   L       : FIR 长度（奇数）
%   method  : 'sinc' | 'lagrange'   （Kaiser-窗sinc 或 Lagrange）
%   beta    : Kaiser β（method='sinc' 时生效）
%   Fs      : 采样率 (Hz)
%   B       : 有效带宽 (Hz)，用于生成带限测试信号与评估
%   N       : 测试信号长度（样点）
%   doplot  : true/false 是否绘图
%
% OUTPUT:
%   stats struct:
%     .mu_set, .mu_int, .mu_frac, .mu_hat, .mu_err, .mu_err_time,
%     .rmse_vs_ideal, .mag_ripple_db([min,max]), .passband([fmin,fmax]/(Fs/2))
%
% NOTE:
%   - 估计器采用“线性互相关 + 抛物线亚样点细化”（稳健、无NaN问题）。
%   - 为与“线性时移”（非循环）对比，提供 ideal_fracdelay_linear()。

if nargin==0
    % ===== 默认示例 =====
    mu      = 0.35;             % 目标分数时延(样点)
    L       = 65;               % FIR 长度（奇数）
    method  = 'sinc';           % 或 'lagrange'
    beta    = 6;                % Kaiser β（~60 dB 旁瓣）
    Fs      = 1e6;              % 采样率
    B       = 2e5;              % 有效带宽
    N       = 65536;            % 样点数
    doplot  = true;             % 绘图开关
end

% ===== 防御与固定参数 =====
assert(mod(L,2)==1, 'L must be odd');
n0 = (L-1)/2;
if B >= Fs, error('B 必须 < Fs'); end

% ===== 生成复带限测试信号 =====
x = bandlimited_noise(N, Fs, B);     % 均值≈0, 方差≈1
x = x(:);

% ===== 目标时延拆分（整数 + 分数） =====
mu_set  = mu;
mu_int  = fix(mu_set);               % 向零取整（也可改 floor/round）
mu_frac = mu_set - mu_int;           % 交给 FIR 的分数残差（建议 |mu_frac|<=0.5）

% ===== 设计 FIR 分数时延核并应用（线性相位，移除群延）=====
h  = design_fd_fir(mu_frac, L, method, beta);
y  = filter(h, 1, [x; zeros(L-1,1)]);
y  = y(n0+1 : n0+N);                 % 去掉 FIR 固有群延 n0

% ===== 再施加整数样点移位（线性、非循环）=====
if     mu_int > 0
    y = [zeros(mu_int,1); y(1:end-mu_int)];
elseif mu_int < 0
    k = -mu_int;
    y = [y(1+k:end); zeros(k,1)];
end

% ===== 稳健估计实际时延（样点）=====
[mu_hat, lags, rxy] = estimate_delay_parabolic(x, y);

% ===== 与“理想线性时移”输出做波形误差对比 =====
y_ideal = ideal_fracdelay_linear(x, mu_set, 1024);  % 大护带近似线性
rmse_vs_ideal = norm(y - y_ideal) / (norm(y_ideal) + eps);

% ===== FIR 幅频（用于看通带起伏）=====
[H, fH] = freqz(h, 1, 4096, Fs);
pass = (abs(fH) <= B/2);
mag_db = 20*log10(abs(H(pass)) + 1e-12);
ripple = [min(mag_db), max(mag_db)];
pb_norm = [min(fH(pass)) max(fH(pass))]/(Fs/2);

% ===== 汇总输出 =====
stats = struct('mu_set',mu_set,'mu_int',mu_int,'mu_frac',mu_frac, ...
               'mu_hat',mu_hat,'mu_err',mu_hat-mu_set, ...
               'mu_err_time',(mu_hat-mu_set)/Fs, ...
               'rmse_vs_ideal',rmse_vs_ideal, ...
               'mag_ripple_db',ripple, ...
               'passband', pb_norm);

% ===== 打印与绘图 =====
fprintf('set mu = %.6f,  hat mu = %.6f,  err = %.6f samp (%.3f us)\n', ...
    stats.mu_set, stats.mu_hat, stats.mu_err, 1e6*stats.mu_err_time);

if doplot
    t = (0:N-1)/Fs;
    figure('Name','FIR Fractional Delay Self-Test','Color','w','Position',[60 60 1200 480]);

    subplot(1,3,1);
    plot(lags, abs(rxy), 'LineWidth', 0.9); grid on;
    xlim([mu_set-10, mu_set+10]); xlabel('Lag (samples)'); ylabel('|R_{xy}|');
    title(sprintf('Cross-corr (peak ~ %.3f samp)', mu_hat));

    subplot(1,3,2);
    plot(fH/1e3, 20*log10(abs(H)+1e-12), 'LineWidth', 1); grid on;
    xline(B/2/1e3,'--r','B/2'); xline(-B/2/1e3,'--r');
    xlabel('Frequency (kHz)'); ylabel('|H| (dB)');
    title(sprintf('%s FIR (L=%d), ripple ~ [%.2f, %.2f] dB', upper(method), L, ripple(1), ripple(2)));

    subplot(1,3,3);
    plot(t, real(y_ideal), 'k-', 'LineWidth', 0.8); hold on;
    plot(t, real(y),       'r-', 'LineWidth', 0.8);
    legend('ideal linear','FIR','Location','best'); grid on;
    xlabel('Time (s)');
    title(sprintf('\\mu_{set}=%.3f, \\mu_{hat}=%.3f, RMSE=%.2e', mu_set, mu_hat, rmse_vs_ideal));
end

end
% ======================= 子函数区 =======================

function h = design_fd_fir(mu, L, method, beta)
% 设计分数时延 FIR：y[n] ≈ x[n - mu] （只实现分数部分，整数部分外部处理）
% L 必须为奇数；mu 可任意实数（常用在 [-0.5,0.5]）
    n0 = (L-1)/2; 
    n  = (0:L-1).';
    switch lower(method)
        case 'sinc'      % Kaiser-窗 sinc
            h = sinc((n - n0) - mu) .* kaiser(L, beta);
            h = h / sum(h);                   % DC 归一（幅度不失真）
        case 'lagrange'  % 拉格朗日（最大平坦群时延）
            k = (0:L-1).';
            h = ones(L,1);
            D = mu + n0;                      % 等效在 [0..L-1] 上的延时
            for i = 1:L
                idx = k; idx(i) = [];
                h(i) = prod( (D - idx) ./ ((i-1) - idx) );
            end
        otherwise
            error('未知 FIR method: %s', method);
    end
    h = double(h(:));
end

function [d_hat, lags, r] = estimate_delay_parabolic(x, y)
% 线性互相关 + 抛物线插值，稳健估计 y≈x 延时（样点；正值≈y滞后x）
    x = double(x(:)); y = double(y(:));
    N = numel(x);
    % 线性互相关（非循环）
    r   = xcorr(x, y);           % 长度 2N-1
    lags= (-N+1):(N-1);
    [~, idx] = max(abs(r));
    if idx>1 && idx<numel(r)
        a = abs(r(idx-1)); b = abs(r(idx)); c = abs(r(idx+1));
        denom = (a - 2*b + c);
        if abs(denom) < eps
            delta = 0;
        else
            delta = 0.5 * (a - c) / denom;    % ∈[-0.5,0.5]
        end
    else
        delta = 0;
    end
    d_hat = lags(idx) + delta;
end

function y = ideal_fracdelay_linear(x, D, G)
% 大护带 + 频域相位斜坡 → 取中段，近似线性（非循环）分数时移
    if nargin<3, G = 1024; end
    x = x(:); N = numel(x);
    L = 2^nextpow2(N + 2*G);
    xx = zeros(L,1,'like',x);
    xx(G+1:G+N) = x;
    k  = (0:L-1).';
    Y  = fft(xx);
    W  = exp(-1j*2*pi*(k/L)*D);
    yy = ifft(Y .* W);
    y  = yy(G+1:G+N);
end

function x = bandlimited_noise(N, Fs, B)
% 生成复带限白噪声，通带 |f|<=B/2
    X = (randn(N,1) + 1j*randn(N,1))/sqrt(2);
    f = ((-N/2):(N/2-1)).' * (Fs/N);
    mask = abs(f) <= (B/2);
    Xf = fftshift(fft(X));
    Xf(~mask) = 0;
    x  = ifft(ifftshift(Xf));
    x  = x / (rms(x) + eps);
end
