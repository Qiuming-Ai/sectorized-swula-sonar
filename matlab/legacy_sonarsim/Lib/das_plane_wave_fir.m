function [img, beam, ranges_out] = das_plane_wave_fir(cfg, rx)
% Plane-wave DAS using only FIR fractional delay (per-angle, no FFT)
% 整数移位 + FIR 分数时延（sinc/lagrange）+（可选）基带载频相位补偿 → 求和
%
% Inputs:
%   cfg.t0, cfg.fs, cfg.c, cfg.rx_xyz, cfg.angles [, cfg.ranges]
%   cfg.fd.sign     : +1/-1（默认 +1；主瓣方向反了可设 -1）
%   cfg.w           : 阵元固定权重 [M×1]（默认全 1）
%   % —— FIR 分数时延参数 —— 
%   cfg.fir.method  : 'sinc' | 'lagrange'   （默认 'sinc'）
%   cfg.fir.L       : FIR 长度（奇数，默认 65）
%   cfg.fir.beta    : Kaiser β（method='sinc' 时生效，默认 6）
%   cfg.fir.guard   : 每侧零护带样点数（默认 max(512, round(0.1*Nsamp))）
%   % —— 基带载频相位补偿（可选）——
%   cfg.f0          : 下变频中心频率（Hz）；提供则视作“基带数据”
%   cfg.mix_sign    : +1 表示下变频用 exp(+j2πf0 t)，-1 表示 exp(-j2πf0 t)（默认 +1）
%
% Outputs:
%   img        : 若提供 cfg.ranges，则为 [Nr×Na] 取样后的角-距图；否则 []
%   beam       : [Lsig×Na] 每角度一整条波束时序（与 ranges_out 匹配）
%   ranges_out : [Lsig×1] 与 beam 第一维匹配的“距离轴”：r = c*(t0 + (0:Lsig-1)/fs)
%
% 说明：
% - 全过程不用 FFT；FIR 线性相位，统一在末尾以 ranges 取样（可选）。
% - 若输入为“基带解调后数据”，强烈建议提供 cfg.f0 与 cfg.mix_sign，以补偿载频相位项。

% ---------- 基本参数 ----------
if ~isfield(cfg,'c')  || isempty(cfg.c),  cfg.c  = 1500;  end
if ~isfield(cfg,'t0') || isempty(cfg.t0), cfg.t0 = 0;     end
t0 = cfg.t0; fs = cfg.fs; c = cfg.c;

angles = -cfg.angles(:).';                      % 与既有实现保持一致：取反
has_ranges = isfield(cfg,'ranges') && ~isempty(cfg.ranges);

[Nsamp, M] = size(rx);
rx_xyz = cfg.rx_xyz;
if size(rx_xyz,2)==2, rx_xyz = [rx_xyz, zeros(M,1)]; end
rx_xyz = double(rx_xyz);
Na = numel(angles);

% ---------- 方向符号 / 权重 ----------
if ~isfield(cfg,'fd') || ~isfield(cfg.fd,'sign') || isempty(cfg.fd.sign)
    sgn = +1;
else
    sgn = sign(cfg.fd.sign);
end
if isfield(cfg,'w') && ~isempty(cfg.w)
    w = cfg.w(:);
    if numel(w) ~= M, error('cfg.w 大小应为 M 元素'); end
else
    w = ones(M,1);
end

% ---------- FIR 参数 ----------
method = getfield_def(cfg,'fir','method','sinc');    % 'sinc' | 'lagrange'
L      = getfield_def(cfg,'fir','L',65);             % 奇数
if mod(L,2)==0, error('cfg.fir.L 必须为奇数'); end
beta   = getfield_def(cfg,'fir','beta',6);           % Kaiser β
G      = getfield_def(cfg,'fir','guard', max(512, round(0.1*Nsamp)) );

% ---------- 基带载频相位补偿开关（可选） ----------
has_f0   = isfield(cfg,'f0') && ~isempty(cfg.f0);
mix_sign = getfield_def(cfg,'','mix_sign', +1);      % +1: 用了 +j2πf0t 下变频；-1: 用了 -j2πf0t

% ================================================================
% ① 全角度预计算样点时延范围 → 统一整数移位不截断的零填充 + 护带
% d_samp_all = sgn * (fs/c) * (x_m*sinθ + z_m*cosθ)
% ================================================================
ux = sin(angles);                 % 1×Na
uz = cos(angles);                 % 1×Na
dx = rx_xyz(:,1);                 % M×1
dz = rx_xyz(:,3);                 % M×1
d_samp_all = sgn * (fs/c) * (dx*ux + dz*uz);  % M×Na

min_int = floor(min(d_samp_all(:)));
max_int = ceil(max(d_samp_all(:)));

pad_left  = max(0, -min_int);
pad_right = max(0,  max_int);

Lsig = Nsamp + pad_left + pad_right + 2*G;    % 统一帧长（移位后不截断 + 左右护带）

% ---------- 输出轴（与 beam 匹配） ----------
ranges_out = c * ( t0 + ((0:Lsig-1).'/fs) )/2;  % 一程距离；若想双程，末尾再 /2

% ---------- 预分配 ----------
beam = complex(zeros(Lsig, Na), 0);

% ================================================================
% ② 逐角度：整数移位 + FIR 分数时延（+可选载频相位）→ 求和
% ================================================================
n0 = (L-1)/2;                                  % FIR 群时延（线性相位）
for ia = 1:Na
    % —— 当前角度的样点域时延 —— 
    d_samp = sgn * (fs/c) * (rx_xyz * [ux(ia); 0; uz(ia)]);  % [M×1]
    n_int  = round(d_samp);                                  % 整数移位（前移/后移）
    mu     = d_samp - n_int;                                 % 分数残差（可在 [-0.5,0.5]）

    % —— 可选：基带载频相位补偿（每通道常相位）——
    if has_f0
        tau   = d_samp / fs;                                 % 秒
        % 若下变频用了 +j2πf0t，残留相位为 +j2πf0 τ → 乘 e^{-j2πf0 τ}
        % 若下变频用了 -j2πf0t，残留相位为 -j2πf0 τ → 乘 e^{+j2πf0 τ}
        rot = exp(-1j * mix_sign * 2*pi * cfg.f0 * tau).';  % 1×M
    else
        rot = ones(1,M,'like',rx);
    end

    % —— 按统一帧长与护带，做整数移位（非循环，零填充）——
    x_shift = zeros(Lsig, M, 'like', rx);
    starts  = G + pad_left + n_int(:) + 1;
    ends    = starts + Nsamp - 1;
    for m = 1:M
        x_shift(starts(m):ends(m), m) = rx(:,m) .* rot(m);   % 先乘常相位再移位
    end

    % —— 每通道 FIR 分数时延（sinc/lagrange），补偿 mu —— 
    % 线性相位 FIR 会引入 n0 样点延迟：用移位消掉保证对齐
    for m = 1:M
        h = design_fd_fir(mu(m), L, method, beta);
        y = filter(h, 1, x_shift(:,m));                     % 长度 = Lsig
        % 去群延（非循环）：把 y 左移 n0
        if n0>0
            x_shift(:,m) = [y(n0+1:end); zeros(n0,1,'like',y)];
        else
            x_shift(:,m) = y;
        end
    end

    % —— 阵元加权求和 → 该角的一条波束时序 —— 
    beam(:, ia) = x_shift * w;
end

% ================================================================
% ③ 若提供 ranges，仅用于最后取样生成 img
% ================================================================
if has_ranges
    ranges = cfg.ranges(:);
    Nr     = numel(ranges);
    idx_r  = round( (ranges./c - t0) * fs ) + 1;            % [Nr×1]
    valid  = (idx_r>=1) & (idx_r<=Lsig);

    img = complex(zeros(Nr, Na));
    ir  = idx_r(valid);
    for ia = 1:Na
        col = complex(zeros(Nr,1),0);
        col(valid) = beam(ir, ia);
        img(:, ia) = col;
    end
else
    img = [];
end
end

% ======================= 辅助函数们 ===========================
function v = getfield_def(s, sub, name, defaultv)
% 读 cfg.sub.name，若无则 defaultv
    if isempty(sub)
        if isfield(s,name) && ~isempty(s.(name))
            v = s.(name);
        else
            v = defaultv;
        end
    else
        if isfield(s,sub) && isfield(s.(sub),name) && ~isempty(s.(sub).(name))
            v = s.(sub).(name);
        else
            v = defaultv;
        end
    end
end

function h = design_fd_fir(mu, L, method, beta)
% 设计分数时延 FIR：实现 y[n] ≈ x[n - mu]
% L 必须为奇数；mu 可为任意实数（常用在 [-0.5,0.5]）
    n0 = (L-1)/2; n = (0:L-1).';
    switch lower(method)
        case 'sinc'      % Kaiser-窗 sinc
            h = sinc((n - n0) - mu) .* kaiser(L, beta);
            h = h / sum(h);                   % DC 归一（幅度不失真）
        case 'lagrange'  % 拉格朗日（最大平坦群时延）
            % 将插值中心放在 n0，构造关于 mu 的 L 点插值核
            k = (0:L-1).';
            h = ones(L,1);
            for i = 1:L
                idx = k; idx(i) = [];         % 0..L-1 但以 n0 为中心对称
                h(i) = prod( ((mu + n0) - idx) ./ ((i-1) - idx) );
            end
        otherwise
            error('未知 cfg.fir.method: %s', method);
    end
end
