%% Peak ratio in dB with WIDEBAND Doppler (time-scaling), ref template at fr (no Doppler)
% Window: Hamming; T=0.5 ms; Bc=20 kHz; variable Δf/Bc; v = [0 1 2 5 10] m/s
% clear; clc;

% -------- Parameters --------
T    = 0.5e-3;                  % pulse width (s)
Bc   = 20e3;                    % subband LFM bandwidth (Hz)
mu   = Bc / T;                  % LFM chirp rate (Hz/s)
c    = 1500;                    % sound speed (m/s)
fr   = 450e3;                   % reference subband center (Hz) <-- set as needed
fs   = 2e6;                     % sample rate (Hz)
t    = (0:1/fs:T-1/fs).';       % time vector [0, T)
N    = numel(t);
M      = 1001;                          % 横轴采样数
vel_list = [0 1 2 5 10];        % m/s (positive = approaching)
df_over_Bc = linspace(0, 2.5, M);
df_vec     = df_over_Bc * Bc;   % Δf in Hz (use +Δf; ±对称)
%% 定义不同时间窗（列向量）
W = struct();
W(1).name = 'Rect';     W(1).w = ones(N,1);
W(2).name = 'Hann';     W(2).w = hann(N);
W(3).name = 'Hamming';  W(3).w = hamming(N);
W(4).name = 'Blackman'; W(4).w = blackman(N);
W(5).name = 'Kaiser\beta=6'; W(5).w = kaiser(N, 6);
% 不同窗的线型/点型
styles  = {'-','--','-.',':','-'};
markers = {'o','s','^','d','v'};
mkEvery = round(M/25);   % 每隔多少点画一个 marker
% -------- Continuous-time Hamming window on [0,T] --------
w_cont = @(tau) ((tau>=0)&(tau<=T)) .* (cos(2*pi*tau/T));

% -------- Reference (no Doppler) baseband template at fr --------
w_t   = w_cont(t);
s_ref = w_t .* exp(1j*pi*mu*(t.^2));  % baseband LFM at fr (no Doppler)
h_mf  = conj(flipud(s_ref));          % matched-filter template

% -------- Wideband Doppler factor (two-way, exact) --------
alpha = @(v) (c + v) ./ (c - v);      % >1 for approaching, <1 for receding

% -------- Doppler-affected baseband echoes (demod at fr) --------
% Passband echo for subband center f_c = fr + df: s_pass(t) ~ s(α t)
% After demod by exp(-j2π fr t), baseband becomes:
%   s_base(t) = w(α t) * exp{ j 2π [ (α-1)fr * t + α*df * t + 0.5*μ*α^2 * t^2 ] }
make_echo_base = @(df, v) ...
    (w_cont(alpha(v)*t)) .* ...
    exp(1j*2*pi*((alpha(v)-1)*fr .* t + alpha(v)*df .* t + 0.5*mu*(alpha(v)^2)*(t.^2)));

% -------- Compute peak ratios (linear) --------
R = zeros(numel(vel_list), numel(df_vec));
for iv = 1:numel(vel_list)
    v   = vel_list(iv);
    s_fr = make_echo_base(0, 0).* W(iv).w;          % fr-band, wideband Doppler
    y_fr = conv(s_fr, h_mf.*W(iv).w, 'full');
    peak_fr = max(abs(y_fr));             % denominator peak (ref template at fr)

    for idf = 1:numel(df_vec)
        df = df_vec(idf);
        s_df = make_echo_base(df, 0).* W(iv).w;     % (fr+Δf)-band, wideband Doppler
        y_df = conv(s_df, h_mf.* W(iv).w, 'full');
        peak_df = max(abs(y_df));
        R(iv, idf) = peak_df / peak_fr;   % linear ratio
    end
end

% -------- Convert to dB (amplitude) --------
dB_floor = -100;                            % floor to avoid -Inf
R_safe   = max(R, 10^(dB_floor/20));        % clip very small values
R_dB     = 20*log10(R_safe);                % amplitude ratio -> dB

% -------- Plot (dB) --------
figure('Color','w'); hold on;
for iv = 1:numel(vel_list)
    plot(df_over_Bc, R_dB(iv,:),styles{iv}, 'LineWidth', 1.6,'Marker', markers{iv}, 'MarkerIndices', 1:mkEvery:M,'DisplayName', W(iv).name);
end
grid on; box on;
ylim([dB_floor 1]);                         % show from floor to ~0 dB
xlabel('\Delta f / B_c');
ylabel('Peak amplitude ratio (dB, ref template at f_r)');
title(sprintf('Wideband Window, T=%.1f ms, B_c=%.0f kHz, f_r=%.0f kHz', ...
      T*1e3, Bc/1e3, fr/1e3));
legend('Location','northeast');
ylim([-80 10]);
% -------- Optional: export --------
Result.df_over_Bc = df_over_Bc;
Result.vel        = vel_list;
Result.R_linear   = R;        % linear ratio
Result.R_dB       = R_dB;     % dB ratio
assignin('base','Result',Result);
