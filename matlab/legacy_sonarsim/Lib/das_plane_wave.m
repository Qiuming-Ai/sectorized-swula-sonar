function img = das_plane_wave(cfg, rx)
% DAS (Delay-and-Sum) for plane-wave transmit, integer-delay only (no frac. interp)
% + SVA (Spatially Variant Apodization, min-processing)
% + CF-family weights: CF / PCF / GCF / MCF / SLSC  [可选像素级乘性权重]
% ------------------------------------------------------------------------
% Inputs:
%   cfg.t0, cfg.fs, cfg.rx_xyz, cfg.angles, cfg.ranges, cfg.c 同原函数
%
%   % —— SVA 相关（同原）——
%   cfg.sva.enable     : true/false（默认 false）
%   cfg.sva.windows    : {'rect','hann','hamming','kaiser:6'} 等
%   cfg.sva.ref_idx    : 相位参考窗索引（默认 1）
%   cfg.sva.norm       : 'mean'（默认）或 'none'
%
%   % —— 相干权重（新增/扩展）——
%   cfg.cf.enable  : true/false（默认 false）
%   cfg.cf.mode    : 'cf' | 'pcf' | 'gcf' | 'mcf' | 'slsc'   % 默认 'cf'
%   cfg.cf.gamma   : >0 权重指数（默认 1.0）
%   cfg.cf.eps     : 稳定项（默认 1e-12）
%   % GCF:
%   cfg.cf.p       : GCF 幂指数 p（默认 2；p=2 退化为 CF），1<p<2温和，p>2要求更高相关性
%   % SLSC:
%   cfg.cf.L       : 短延迟最大滞后（默认 min(16, floor(M/4))）
%
%   rx             : 原始信号 [Nsamp×M]，复/实均可
%
% Output:
%   img            : 复值图 [Nr×Na]
%
% ---------------- 参数与输入检查 ----------------
if ~isfield(cfg,'c') || isempty(cfg.c), cfg.c = 1500; end
t0     = cfg.t0;
fs     = cfg.fs;
c      = cfg.c;
angles = -cfg.angles(:).';     % 与原函数保持一致
ranges = cfg.ranges(:);

[Nsamp, M] = size(rx);

rx_xyz = cfg.rx_xyz;
if size(rx_xyz,2) == 2
    rx_xyz = [rx_xyz, zeros(M,1)]; % 补 y=0
end
rx_xyz = double(rx_xyz);

Na = numel(angles);
Nr = numel(ranges);

% ---------------- SVA 窗集合准备 ----------------
if ~isfield(cfg,'sva') || ~isfield(cfg.sva,'enable') || ~cfg.sva.enable
    sva_enable = false;
    W = ones(M,1); Kwin = 1; ref_idx = 1; normMode = 'mean';
else
    sva_enable = true;
    if ~isfield(cfg.sva,'windows') || isempty(cfg.sva.windows)
        winSpec = {'rect','hann','hamming'};
    else
        winSpec = cfg.sva.windows;
    end
    if ~isfield(cfg.sva,'ref_idx') || isempty(cfg.sva.ref_idx)
        ref_idx = 1;
    else
        ref_idx = max(1, min(numel(winSpec), round(cfg.sva.ref_idx)));
    end
    if ~isfield(cfg.sva,'norm') || isempty(cfg.sva.norm)
        normMode = 'mean';
    else
        normMode = cfg.sva.norm;
    end
    W = build_windows(M, winSpec, normMode);    % [M×K]
    Kwin = size(W,2);
end

% ---------------- CF 家族选项 ----------------
if ~isfield(cfg,'cf') || ~isfield(cfg.cf,'enable') || ~cfg.cf.enable
    cf_enable = false;
    cf_mode   = 'cf';
    cf_gamma  = 1.0;
    cf_eps    = 1e-12;
    gcf_p     = 2;
    slsc_L    = 0;
else
    cf_enable = true;
    cf_mode   = getfield_default(cfg.cf,'mode','cf');   % 'cf' | 'pcf' | 'gcf' | 'mcf' | 'slsc'
    cf_gamma  = getfield_default(cfg.cf,'gamma',1.0);
    cf_eps    = getfield_default(cfg.cf,'eps',1e-12);
    gcf_p     = getfield_default(cfg.cf,'p',2);
    slsc_L    = getfield_default(cfg.cf,'L', min(16, max(1,floor(M/4))));
    slsc_L    = min(slsc_L, max(0,M-1));  % L ≤ M-1
end

% 依据模式准备额外统计
need_pcf  = cf_enable && strcmpi(cf_mode,'pcf');
need_gcf  = cf_enable && strcmpi(cf_mode,'gcf');
need_mcf  = cf_enable && strcmpi(cf_mode,'mcf');
need_cf   = cf_enable && strcmpi(cf_mode,'cf');
need_slsc = cf_enable && strcmpi(cf_mode,'slsc') && slsc_L>0;

img = complex(zeros(Nr, Na));   % 输出

% ---------------- 主循环（按角度） ----------------
for ia = 1:Na
    alpha = angles(ia);
    sa = sin(alpha);  ca = cos(alpha);

    % 本角度下，所有 range 的像素坐标 [Nr×3]
    px = ranges .* sa;             % x
    py = zeros(Nr,1);              % y（2D成像）
    pz = ranges .* ca;             % z
    tau_tx = ranges ./ c;          % 平面波 TX 延时 [Nr×1]

    % 计算到各通道的 RX 延时矩阵 tau_rx [Nr×M]
    dx = px - rx_xyz(:,1).';       % [Nr×M]
    dy = py - rx_xyz(:,2).';
    dz = pz - rx_xyz(:,3).';
    R  = sqrt(dx.^2 + dy.^2 + dz.^2);   % [Nr×M]
    tau_rx = R ./ c;                    % [Nr×M]

    % 总延时 τ_total = τ_tx + τ_rx
    tau_total = tau_rx + tau_tx;        % [Nr×M]
    idx   = round( (tau_total - t0) * fs ) + 1;   % [Nr×M]，1-based
    valid = (idx >= 1) & (idx <= Nsamp);

    % ---- 多窗累加（SVA） ----
    accSVA   = complex(zeros(Nr, Kwin));  % [Nr×K]

    % ---- 通用统计（CF/GCF/MCF/PCF） ----
    if cf_enable
        s_coh    = complex(zeros(Nr,1));      % Σ x_m
        s_abs2   = zeros(Nr,1);               % Σ |x_m|^2
        s_abs1   = zeros(Nr,1);               % Σ |x_m|   (MCF)
        if need_gcf
            s_absp = zeros(Nr,1);             % Σ |x_m|^p (GCF)
        end
        cnt_valid= zeros(Nr,1);               % 有效通道数
        if need_pcf
            ph_sum = complex(zeros(Nr,1));    % Σ e^{j∠x_m}
        end
    end

    % ---- SLSC 统计（短延迟归一化互相关的幅度平均）----
    if need_slsc
        slsc_acc = zeros(Nr, slsc_L);         % 每个滞后的累计相关
        slsc_cnt = zeros(Nr, slsc_L);         % 每个滞后的有效对数
        buf      = complex(zeros(Nr, slsc_L));  % 环形缓冲：存前 L 个通道
        buf_abs  = zeros(Nr, slsc_L);
        buf_val  = false(Nr, slsc_L);
        buf_pos  = 0;                          % 写入指针（1..L）
    end

    % 逐通道取样并累加
    for m = 1:M
        im = idx(:,m);
        vmask = valid(:,m);
        vm = complex(zeros(Nr,1,'like',rx));
        if any(vmask)
            imv = im(vmask);
            vm(vmask) = rx(imv, m); % 该通道对所有 range 的取样
        end

        % —— SVA 累加： accSVA += vm * W(m,:) —— 
        accSVA = accSVA + bsxfun(@times, vm, W(m,:));

        % —— 通用 CF-family 统计 —— 
        if cf_enable
            s_coh   = s_coh   + vm;           % Σ x_m
            s_abs2  = s_abs2  + abs(vm).^2;   % Σ |x_m|^2
            s_abs1  = s_abs1  + abs(vm);      % Σ |x_m|
            cnt_valid = cnt_valid + double(vmask);
            if need_gcf
                s_absp = s_absp + abs(vm).^gcf_p;
            end
            if need_pcf
                tmp = complex(zeros(Nr,1));
                idxv = find(vmask);
                if ~isempty(idxv)
                    xv = vm(idxv);
                    tmp(idxv) = exp(1j*angle(xv));
                end
                ph_sum = ph_sum + tmp;
            end
        end

        % —— SLSC: 与前面 L 个通道做短延迟相关（流式，不存 NxM）——
        if need_slsc
            buf_pos = mod(buf_pos, slsc_L) + 1;  % 写入当前通道的位置
            % 计算当前 vm 与 m-k（k=1..L）的归一化互相关幅度
            if any(vmask)
                abs_vm = abs(vm);
                for k = 1:slsc_L
                    idxk = buf_pos - k;  if idxk<=0, idxk = idxk + slsc_L; end
                    pair_mask = vmask & buf_val(:, idxk);
                    if any(pair_mask)
                        prev   = buf(pair_mask, idxk);
                        prevA  = buf_abs(pair_mask, idxk);
                        denom  = abs_vm(pair_mask) .* prevA + 1e-12;
                        cohmag = abs( vm(pair_mask) .* conj(prev) ) ./ denom; % ∈[0,1]
                        slsc_acc(pair_mask, k) = slsc_acc(pair_mask, k) + cohmag;
                        slsc_cnt(pair_mask, k) = slsc_cnt(pair_mask, k) + 1;
                    end
                end
            end
            % 更新环形缓冲
            buf(:,      buf_pos) = vm;
            buf_abs(:,  buf_pos) = abs(vm);
            buf_val(:,  buf_pos) = vmask;
        end
    end

    % —— 形成 SVA 或常规 DAS 输出（复值）——
    if sva_enable
        y_ref   = accSVA(:, ref_idx);      % 参考窗相位
        mags    = abs(accSVA);             % [Nr×K]
        minMag  = min(mags, [], 2);        % [Nr×1]
        phase_r = exp(1j * angle(y_ref));  % [Nr×1]
        y_pix   = phase_r .* minMag;       % [Nr×1]
    else
        y_pix   = accSVA(:,1);             % 矩形窗 DAS
    end

    % —— 相干权重（像素级乘性）——
    if cf_enable
        Meff = max(cnt_valid, 1);  % 每像素有效通道数
        switch lower(cf_mode)
            case 'cf'   % |Σx|^2 / (Meff * Σ|x|^2)
                cf_fac = (abs(s_coh).^2) ./ (Meff .* s_abs2 + cf_eps);

            case 'pcf'  % | (1/Meff) * Σ e^{j∠x} |
                cf_fac = abs(ph_sum) ./ Meff;

            case 'gcf'  % (|Σx|^p) / (Meff^(p-1) * Σ|x|^p)
                num = abs(s_coh).^gcf_p;
                den = (Meff.^(max(gcf_p-1,0))) .* (s_absp + cf_eps);
                cf_fac = num ./ den;

            case 'mcf'  % |Σx|^2 / ( (Σ|x|)^2 )
                den = (s_abs1.^2) + cf_eps;
                cf_fac = (abs(s_coh).^2) ./ den;

            case 'slsc' % mean_k 〈 | x_m x_{m-k}* | / (|x_m||x_{m-k}|) 〉
                % 平均每个滞后的相关，再对 k=1..L 做平均
                slsc_mean = zeros(Nr,1);
                if slsc_L>0
                    % 防除零
                    slsc_cnt = max(slsc_cnt, 1);
                    Ck = slsc_acc ./ slsc_cnt;          % [Nr×L]，每个滞后的相关幅
                    slsc_mean = mean(Ck, 2);            % [Nr×1]
                end
                cf_fac = slsc_mean;

            otherwise
                error('未知 cfg.cf.mode: %s', cf_mode);
        end

        if cf_gamma ~= 1
            cf_fac = cf_fac .^ cf_gamma;
        end
        % 裁剪到 [0,1]，避免数值溢出（gcf/mcf 在边界处可能略>1）
        cf_fac = min(max(real(cf_fac), 0), 1);

        y_pix = y_pix .* cf_fac;
    end

    img(:, ia) = y_pix;
end
end

% ================== 辅助：生成窗集合 ==================
function W = build_windows(M, winSpec, normMode)
    if isnumeric(winSpec)
        W = double(winSpec);
        if size(W,1) ~= M
            error('数值窗矩阵行数必须等于 M');
        end
    else
        if ischar(winSpec), winSpec = {winSpec}; end
        K = numel(winSpec);
        W = zeros(M, K);
        n = (0:M-1).';
        for k = 1:K
            spec = winSpec{k};
            if ischar(spec)
                if strncmpi(spec,'kaiser',6)
                    beta = 6;
                    tok = strsplit(spec,':');
                    if numel(tok)>=2
                        beta = str2double(tok{2});
                        if ~isfinite(beta), beta = 6; end
                    end
                    x = (2*n)/(M-1) - 1;
                    w = besseli(0, beta*sqrt(max(0,1 - x.^2))) / besseli(0, beta);
                elseif strcmpi(spec,'rect')
                    w = ones(M,1);
                elseif any(strcmpi(spec, {'hann','hanning'}))
                    if M==1, w = 1; else
                        w = 0.5 - 0.5*cos(2*pi*n/(M-1));
                    end
                elseif strcmpi(spec,'hamming')
                    if M==1, w = 1; else
                        w = 0.54 - 0.46*cos(2*pi*n/(M-1));
                    end
                else
                    error('未知窗类型: %s', spec);
                end
            else
                error('winSpec 中存在非字符串项');
            end
            switch lower(normMode)
                case 'mean'
                    mu = mean(w);
                    if mu>0, w = w / mu; end   % DC 增益与矩形窗对齐
                case 'none'
                otherwise
                    error('未知归一化模式: %s', normMode);
            end
            W(:,k) = w(:);
        end
    end
end

% ================== 辅助：安全取字段默认值 ==================
function v = getfield_default(s, name, defaultv)
    if isfield(s,name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = defaultv;
    end
end
