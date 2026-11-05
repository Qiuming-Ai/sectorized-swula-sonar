function [y, t0_ref] = sim_rx_from_scatterers_simple(P, A, TX, RX, f0, c, fs, excitation, opts)
% sim_rx_from_scatterers_simple
% 纯整数时延、最简单成分的超声回波合成（Born 近似，单次散射，双程1/R扩散）
%
% 输入
%   P   : [K x 3]  散射体坐标 (m)
%   A   : [K x 1]  散射体幅值/散射系数（线性幅度）
%   TX  : [Nt x 3] 发射阵元坐标 (m)（假定同时同相激励）
%   RX  : [Nr x 3] 接收阵元坐标 (m)
%   f0  : 中心频率 (Hz)（此简化模型不直接使用；保留作记录）
%   c   : 声速 (m/s)
%   fs  : 采样率 (Hz)
%   excitation : [Nx1] 发射信号（列向量）
%   opts.round : 'round'(默认) | 'floor' | 'ceil'  时延→样点的取整方式
%   opts.precision : 'single'(默认) | 'double'
%
% 输出
%   y       : [M x Nr] 各通道接收信号（已与 excitation 卷积；最早到达对齐到 t=0）
%   t0_ref  : 最早到达的绝对时延(秒)，供需要时标的场合使用
%
% 说明
% - 不含吸收/频散/指向性/机电脉冲响应；若需要可在 A 或 excitation 中体现。
% - 复杂度 O(K*Nt*Nr)；K或阵元多时请注意计算量。

    if nargin < 9 || ~isfield(opts,'round'),      opts.round = 'round'; end
    if ~isfield(opts,'precision'),                opts.precision = 'single'; end
    if ~isfield(opts,'delay'),                opts.delay = zeros(size(A)); end

    P  = double(P);   A  = double(A(:));
    TX = double(TX);  RX = double(RX);
    fs = double(fs)*10;  c  = double(c);
    excitation = resample(excitation,10,1);
    K  = size(P,1);
    Nt = size(TX,1);
    Nr = size(RX,1);
    x  = excitation(:);
    Nx = numel(x);

    % ------------------ Pass 1：仅统计最早/最晚到达(确定长度) ------------------
    t_min = inf;  t_max = -inf;
    for r = 1:Nr
        Rr = RX(r,:);
        for k = 1:K
            Pk = P(k,:);
            Rrx = norm(Pk - Rr);                  % 散点→该接收阵元
            % 向量化计算所有 TX→散点距离
            Rtx = sqrt(sum((TX - Pk).^2, 2));    % [Nt x 1]
            tau_vec = (Rtx + Rrx) / c;           % [Nt x 1]
            t_min = min(t_min, min(tau_vec));
            t_max = max(t_max, max(tau_vec)+max(opts.delay));
        end
    end
    t0_ref = t_min;                               % 以最早到达作为时间零点

    % 取整函数
    switch lower(opts.round)
        case 'round', idxfun = @(t) round(t*fs);
        case 'floor', idxfun = @(t) floor(t*fs);
        case 'ceil',  idxfun = @(t) ceil(t*fs);
        otherwise, error('opts.round 仅支持 round/floor/ceil');
    end

    % 冲激响应长度（最晚相对时延处的采样点 + 1）
    Lir = idxfun(t_max - t0_ref) + 1;
    ir  = zeros(Lir, Nr, opts.precision);

    % ------------------ Pass 2：真正累加稀疏冲激响应 ------------------
    for r = 1:Nr
        Rr = RX(r,:);
        for k = 1:K
            Pk  = P(k,:);
            Rrx = norm(Pk - Rr);                  % 散点→接收
            Rtx = sqrt(sum((TX - Pk).^2, 2));     % 所有发射→散点
            tau = (Rtx + Rrx) / c+opts.delay(k);                % [Nt x 1] 绝对时延
            n   = idxfun(tau - t0_ref);           % 相对零点后的采样索引（整数）
            n(n < 0) = 0;                         % 数值保险
            % w   = (A(k) ./ (Rtx .* Rrx));         % 幅度：双程 1/R 衰减
            % w   = sqrt(A(k) ./ (Rtx .* Rrx));
            w   = A(k);
            % 累加到该通道的 IR（同一采样点可能多条路径重合）
            for j = 1:Nt
                ir(n(j)+1, r) = ir(n(j)+1, r) + cast(w(j), opts.precision);
            end
        end
    end

    % ------------------ 与激励卷积，得到通道时域信号 ------------------
    M = Lir + Nx - 1;
    y = zeros(M, Nr, opts.precision);
    x = cast(x, opts.precision);
    for r = 1:Nr
        if any(ir(:,r))
            y(:,r) = conv(x, ir(:,r), 'full');
        end
    end
    y = resample(y,1,10);
end
