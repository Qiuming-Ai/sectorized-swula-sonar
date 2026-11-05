function [BP_dB, theta_deg, Y_all] = mf_beamscan_peakdir(X, H, RX, fs, c, opts)
% mf_beamscan_peakdir
% 多通道 → 逐子带匹配滤波 → 远场整数延时DAS扫描 → 指向性曲线（峰值度量）
%
% 输入
%   X   : [Nt x M]    多通道时间序列
%   H   : [L  x Na]   Na 路匹配滤波器（每列一条，已 conj(flipud(.)) 好）
%   RX  : [Nt x 3]    接收阵元坐标 (m)
%   fs  : 标量        采样率 (Hz)
%   c   : 标量        声速 (m/s)
%   opts.round     : 'round'(默认) | 'floor' | 'ceil'  扫描延时取整
%   opts.theta_deg : [1 x Nθ] 扫描角列表（默认 -85:85）
%   opts.apod      : [Nt x 1] 通道加权（默认全1；也可 'hamming' 自动生成）
%   opts.envelope  : true/false 匹配滤波后是否先取包络做扫描（默认 false）
%   opts.normalize : 'per_subband'(默认) | 'global'  dB 归一方式
%
% 输出
%   BP_dB     : [Na x Nθ]  每子带指向性（dB，0 dB 为各自或全局峰值）
%   theta_deg : [1 x Nθ]   扫描角列表（度）
%   Y_all     : struct     仅供调试/可视化（匹配滤波后数据）
%               .Y{a}      [Nt x My] 第 a 子带匹配滤波后的通道数据
%
% 说明
% - 远场模型：s(θ) = [sinθ, 0, cosθ]，τ_j(θ) = dot(RX(j,:), s)/c
% - 扫描度量：对齐求和后的 |y_sum(t)| 的全时域峰值（max over t）
% - 复杂度：O(Na * (Nt*conv + Nt*Nθ*对齐求和))

    arguments
        X double
        H double
        RX double
        fs double {mustBePositive}
        c  double {mustBePositive}
        opts.round char {mustBeMember(opts.round,{'round','floor','ceil'})} = 'round'
        opts.theta_deg double = -85:1:85
        opts.apod = []
        opts.envelope logical = false
        opts.normalize char {mustBeMember(opts.normalize,{'per_subband','global'})} = 'per_subband'
    end

    [Nt, M] = size(X);
    [Nt_rx, dim] = size(RX);
    assert(Nt_rx==Nt && dim==3, 'RX 必须是 [Nt x 3]');

    [L, Na] = size(H);
    theta_deg = opts.theta_deg(:).';
    Ntheta = numel(theta_deg);

    % --- 通道加权（apodization） ---
    if isempty(opts.apod)
        w = ones(Nt,1);
    elseif ischar(opts.apod) || (isstring(opts.apod) && isscalar(opts.apod))
        sw = string(opts.apod);
        switch sw
            case "hamming"
                w = hamming(Nt);
            otherwise
                error('未知 apod 选项：%s', char(sw));
        end
    else
        w = opts.apod(:);
        assert(numel(w)==Nt, 'opts.apod 长度必须为 Nt');
    end
    w = w / sum(w); % 归一，避免受阵元个数影响

    % --- 用于整数样点取整 ---
    switch opts.round
        case 'round', idxfun = @(t) round(t*fs);
        case 'floor', idxfun = @(t) floor(t*fs);
        case 'ceil',  idxfun = @(t) ceil(t*fs);
    end

    % 预存输出
    BP = zeros(Na, Ntheta);

    % 匹配滤波器翻转共轭（如果传入是参考信号，可在外部先处理；这里假定已是匹配滤波器）
    % h_a = H(:,a) 直接用于 conv

    % --- 逐子带处理 ---
    Y_all = struct('Y', cell(Na,1));
    for a = 1:Na
        h = H(:,a);                        % [L x 1]
        % 对每个通道做匹配滤波（全卷积）
        % 使用conv2按列卷积：把 X 转置成 [M x Nt]，与 h 做 conv，回转
        Y = conv2(X.', h, 'full').';       % [Nt x (M+L-1)]
        if opts.envelope
            % 可选：先取希尔伯特包络（复解析信号幅值）
            for j = 1:Nt
                Y(j,:) = abs(hilbert(Y(j,:)));
            end
        end
        My = size(Y,2);
        Y_all(a).Y = Y;

        % --- 扫描角延时表（样点），统一移成 >=0 ---
        th = deg2rad(theta_deg);
        S  = [sin(th); zeros(1,Ntheta); cos(th)];   % [3 x Nθ]
        tau = (RX * S) ./ c;                        % [Nt x Nθ]（秒）
        dSamp = idxfun(tau);                        % [Nt x Nθ]（样点，可正负）
        dmin = min(dSamp, [], 'all');
        if dmin < 0
            dSamp = dSamp - dmin;                   % 平移使所有延时非负
        end
        maxDelay = max(dSamp, [], 'all');

        % 逐角度做整数对齐求和，并取峰值作为指向性度量
        for it = 1:Ntheta
            d = dSamp(:,it);                        % [Nt x 1]
            % 将各通道按 d(j) 右移后加权求和
            Lsum = My + maxDelay;                   % 足够长
            ysum = zeros(1, Lsum);
            for j = 1:Nt
                idx = (1:My) + d(j);
                ysum(idx) = ysum(idx) + w(j) * Y(j,:);
            end
            % 峰值度量
            BP(a, it) = max(abs(ysum));
        end
    end

    % --- 归一化并转 dB ---
    if strcmp(opts.normalize,'per_subband')
        % 每个子带各自归一（0 dB 在该子带峰值）
        for a = 1:Na
            m = max(BP(a,:)) + eps;
            BP_dB(a,:) = 20*log10(BP(a,:)/m + eps);
        end
    else
        % 全局归一（0 dB 在所有子带/角度的最大值）
        mg = max(BP,[],'all') + eps;
        BP_dB = 20*log10(BP/mg + eps);
    end
end
