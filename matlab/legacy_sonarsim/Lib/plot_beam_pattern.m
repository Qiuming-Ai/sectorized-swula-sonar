function [beam_dB,beam_response, angles] = plot_beam_pattern(rx_xyz, frequency, varargin)
% PLOT_BEAM_PATTERN 计算并绘制阵列波束图
%
% 输入参数:
%   rx_xyz    - 阵列元素坐标 [N×3] 或 [N×2] (x, y, z) 或 (x, z)
%   frequency - 频率 (Hz)
%
% 可选参数 (通过 varargin 或结构体):
%   'steer_angle'  - 扫描角度 (度或弧度), 默认 0
%   'angle_range'  - 观察角度范围 [min_angle, max_angle] (度), 默认 [-90, 90]
%   'angle_step'   - 角度步长 (度), 默认 0.1
%   'c'            - 声速 (m/s), 默认 1500
%   'weights'      - 阵元权重 [N×1], 默认全1 (均匀加权)
%   'window'       - 窗函数类型 ('rect', 'hann', 'hamming', 'kaiser'), 默认 'rect'
%   'plot_type'    - 绘图类型 ('polar', 'rect', 'both'), 默认 'both'
%   'dyn_range'    - 动态范围 (dB), 默认 60
%   'units'        - 角度单位 ('deg' 或 'rad'), 默认 'deg'
%
% 输出:
%   beam_response - 波束响应 [N_angles×1], 复值
%   angles        - 观察角度 [N_angles×1]
%
% 示例:
%   % 线性阵列
%   rx_xyz = [linspace(-0.1, 0.1, 32)', zeros(32, 2)];
%   plot_beam_pattern(rx_xyz, 450e3, 'steer_angle', 0, 'angle_range', [-60, 60]);
%
%   % 2D阵列
%   [x, y] = meshgrid(linspace(-0.05, 0.05, 8));
%   rx_xyz = [x(:), y(:), zeros(64, 1)];
%   plot_beam_pattern(rx_xyz, 450e3);

    % ========== 参数解析 ==========
    p = inputParser;
    addParameter(p, 'steer_angle', 0, @isnumeric);
    addParameter(p, 'angle_range', [-90, 90], @(x) isnumeric(x) && numel(x)==2);
    addParameter(p, 'angle_step', 0.1, @isnumeric);
    addParameter(p, 'c', 1500, @isnumeric);
    addParameter(p, 'weights', [], @isnumeric);
    addParameter(p, 'window', 'rect', @(x) ischar(x) || isstring(x));
    addParameter(p, 'plot_type', 'both', @(x) ismember(x, {'polar', 'rect', 'both','None'}));
    addParameter(p, 'dyn_range', 60, @isnumeric);
    addParameter(p, 'units', 'deg', @(x) ismember(x, {'deg', 'rad'}));
    parse(p, varargin{:});
    
    steer_angle = p.Results.steer_angle;
    angle_range = p.Results.angle_range;
    angle_step = p.Results.angle_step;
    c = p.Results.c;
    weights = p.Results.weights;
    window_type = p.Results.window;
    plot_type = p.Results.plot_type;
    dyn_range = p.Results.dyn_range;
    units = p.Results.units;
    
    % ========== 输入检查 ==========
    if size(rx_xyz, 2) == 2
        rx_xyz = [rx_xyz, zeros(size(rx_xyz, 1), 1)]; % 补 y=0
    elseif size(rx_xyz, 2) ~= 3
        error('rx_xyz 必须是 [N×2] 或 [N×3] 矩阵');
    end
    
    M = size(rx_xyz, 1); % 阵元数量
    rx_xyz = double(rx_xyz);
    
    % ========== 角度处理 ==========
    if strcmp(units, 'deg')
        % 输入为度，转换为弧度用于计算
        if max(abs(angle_range)) > 2*pi
            % 已经是度
            angles_deg = angle_range(1):angle_step:angle_range(2);
            angles_rad = deg2rad(angles_deg);
            steer_angle_rad = deg2rad(steer_angle);
        else
            % 可能是弧度，但假设是度
            angles_deg = angle_range(1):angle_step:angle_range(2);
            angles_rad = deg2rad(angles_deg);
            steer_angle_rad = deg2rad(steer_angle);
        end
    else
        % 输入为弧度
        angles_rad = angle_range(1):angle_step:angle_range(2);
        angles_deg = rad2deg(angles_rad);
        steer_angle_rad = steer_angle;
    end
    
    angles = angles_deg(:); % 输出用度
    Na = numel(angles_rad);
    
    % ========== 权重处理 ==========
    if isempty(weights)
        weights = ones(M, 1);
    else
        if numel(weights) ~= M
            error('权重向量长度必须等于阵元数量');
        end
        weights = weights(:);
    end
    
    % 应用窗函数
    if ~strcmpi(window_type, 'rect')
        win = get_window(M, window_type);
        weights = weights .* win;
    end
    
    % 归一化权重
    weights = weights / sum(abs(weights));
    
    % ========== 计算波束图 ==========
    % 波数
    k = 2 * pi * frequency / c;
    
    % 扫描方向向量 (单位向量)
    u_steer = [sin(steer_angle_rad); 0; cos(steer_angle_rad)]; % 2D平面: x-z平面
    
    % 计算阵列流形向量和波束响应
    beam_response = complex(zeros(Na, 1));
    
    for ia = 1:Na
        % 观察方向向量 (单位向量)
        u_obs = [sin(angles_rad(ia)); 0; cos(angles_rad(ia))];
        
        % 阵列流形向量: a(θ) = exp(-j*k*r·u)
        % 对于扫描波束: a(θ) = exp(-j*k*r·(u_obs - u_steer))
        phase_diff = k * rx_xyz * (u_obs - u_steer); % [M×1]
        a = exp(-1j * phase_diff); % [M×1]
        
        % 波束响应: B(θ) = w^H * a(θ)
        beam_response(ia) = weights' * a;
    end
    
    % ========== 绘图 ==========
    % 归一化并转换为dB
    beam_power = abs(beam_response).^2;
    beam_power_norm = beam_power / max(beam_power);
    beam_dB = 10 * log10(beam_power_norm + eps);
    beam_dB = max(beam_dB, -dyn_range); % 限制动态范围
    
    if strcmp(plot_type, 'polar') || strcmp(plot_type, 'both')
        figure;
        % 兼容性检查：使用 polarplot (R2016a+) 或 polar (旧版本)
        if exist('polarplot', 'builtin') || exist('polarplot', 'file')
            polarplot(deg2rad(angles), beam_power_norm, 'LineWidth', 1.5);
            rlim([0, 1]);
        else
            % 旧版本使用 polar 函数
            polar(deg2rad(angles), beam_power_norm, 'LineWidth', 1.5);
        end
        title(sprintf('Beam Pattern (f=%.1f kHz, steer=%.1f°)', ...
            frequency/1e3, rad2deg(steer_angle_rad)));
        grid on;
    end
    
    if strcmp(plot_type, 'rect') || strcmp(plot_type, 'both')
        if strcmp(plot_type, 'both')
            figure;
        else
            figure;
        end
        plot(angles, beam_dB, 'LineWidth', 1.5);
        xlabel('Angle (degrees)');
        ylabel('Beam Response (dB)');
        title(sprintf('Beam Pattern (f=%.1f kHz, steer=%.1f°)', ...
            frequency/1e3, rad2deg(steer_angle_rad)));
        grid on;
        ylim([-dyn_range, 5]);
        xlim([angles(1), angles(end)]);
        hold on;
        % 标记扫描方向
        plot([rad2deg(steer_angle_rad), rad2deg(steer_angle_rad)], ...
            [-dyn_range, 5], 'r--', 'LineWidth', 1, 'DisplayName', 'Steer Angle');
        legend('Beam Pattern', 'Steer Angle', 'Location', 'best');
        hold off;
    else

    end
    
end

% ========== 辅助函数：生成窗函数 ==========
function w = get_window(M, window_type)
    n = (0:M-1)';
    switch lower(window_type)
        case 'rect'
            w = ones(M, 1);
        case {'hann', 'hanning'}
            if M == 1
                w = 1;
            else
                w = 0.5 - 0.5 * cos(2*pi*n/(M-1));
            end
        case 'hamming'
            if M == 1
                w = 1;
            else
                w = 0.54 - 0.46 * cos(2*pi*n/(M-1));
            end
        case 'kaiser'
            beta = 6; % 默认beta
            if M == 1
                w = 1;
            else
                x = (2*n)/(M-1) - 1;
                w = besseli(0, beta*sqrt(max(0, 1 - x.^2))) / besseli(0, beta);
            end
        otherwise
            warning('未知的窗函数类型，使用矩形窗');
            w = ones(M, 1);
    end
end

