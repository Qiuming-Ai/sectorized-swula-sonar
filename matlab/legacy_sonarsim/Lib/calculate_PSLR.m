function PSLR = calculate_PSLR(beam_pattern, steer_angle, angle_range, varargin)
% CALCULATE_PSLR 计算波束图的峰值旁瓣比（Peak Sidelobe Level Ratio）
%
% 输入参数:
%   beam_pattern - 波束图数据 [N×1] (dB值)
%   steer_angle  - 扫描角度 (度)，通常为主瓣指向角度
%   angle_range  - 角度范围 [min_angle, max_angle] (度)
%                  或角度向量 [N×1] (度)
%
% 可选参数:
%   'first_null' - 第一零点位置 (度)，如果提供则直接使用
%                  如果不提供，将自动调用 find_first_null 函数计算
%   'null_tolerance' - 第一零点容差 (度)，用于定义主瓣范围，默认使用第一零点
%
% 输出参数:
%   PSLR         - 峰值旁瓣比 (dB)，即 主瓣峰值 - 最高旁瓣峰值
%
% 说明:
%   PSLR (Peak Sidelobe Level Ratio) 定义为：
%   PSLR = 主瓣峰值(dB) - 最高旁瓣峰值(dB)
%   
%   主瓣范围：由第一零点定义，通常为 [-first_null, +first_null]
%   旁瓣范围：主瓣范围外的所有角度
%
% 示例:
%   % 使用角度范围，自动计算第一零点
%   beam_db = ...;  % 波束图数据
%   PSLR = calculate_PSLR(beam_db, 0, [-85, 85]);
%
%   % 使用角度向量，自动计算第一零点
%   angles = -85:0.1:85;
%   beam_db = ...;  % 波束图数据
%   PSLR = calculate_PSLR(beam_db, 0, angles);
%
%   % 指定第一零点位置
%   PSLR = calculate_PSLR(beam_db, 0, angles, 'first_null', 1.8);

    % ========== 参数解析 ==========
    p = inputParser;
    addParameter(p, 'first_null', [], @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'null_tolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, varargin{:});
    
    first_null_user = p.Results.first_null;
    null_tolerance = p.Results.null_tolerance;
    
    % ========== 输入验证 ==========
    beam_pattern = beam_pattern(:);  % 确保是列向量
    N = length(beam_pattern);
    
    % 判断angle_range是范围还是角度向量
    if length(angle_range) == 2 && numel(angle_range) == 2
        % 是角度范围 [min, max]
        min_angle = angle_range(1);
        max_angle = angle_range(2);
        angles = linspace(min_angle, max_angle, N)';
    else
        % 是角度向量
        angles = angle_range(:);  % 确保是列向量
        if length(angles) ~= N
            error('角度向量长度必须与波束图数据长度相同');
        end
    end
    
    % ========== 找到主瓣峰值位置 ==========
    % 在主瓣附近搜索峰值（通常峰值在steer角度附近）
    search_range = 10;  % 在steer角度±10度范围内搜索峰值
    peak_search_mask = abs(angles - steer_angle) <= search_range;
    peak_search_indices = find(peak_search_mask);
    
    if isempty(peak_search_indices)
        % 如果搜索范围太小，在整个范围内搜索
        [~, peak_idx] = max(beam_pattern);
    else
        [~, local_max_idx] = max(beam_pattern(peak_search_indices));
        peak_idx = peak_search_indices(local_max_idx);
    end
    
    main_peak_db = beam_pattern(peak_idx);
    
    % ========== 确定第一零点位置 ==========
    if isempty(first_null_user)
        % 自动计算第一零点
        try
            [left_null, right_null] = find_first_null(beam_pattern, steer_angle, angle_range);
            
            % 计算从steer角度到左右零点的距离，取较大者作为第一零点
            if ~isnan(left_null) && ~isnan(right_null)
                % 计算从steer角度到零点的距离
                dist_left = abs(left_null - steer_angle);
                dist_right = abs(right_null - steer_angle);
                first_null = max(dist_left, dist_right);
            elseif ~isnan(left_null)
                first_null = abs(left_null - steer_angle);
            elseif ~isnan(right_null)
                first_null = abs(right_null - steer_angle);
            else
                % 如果找不到第一零点，使用默认值（基于steer角度）
                warning('无法自动找到第一零点，使用默认估计值 2 度');
                first_null = 2.0;
            end
        catch
            warning('自动计算第一零点失败，使用默认估计值 2 度');
            first_null = 2.0;
        end
    else
        % 使用用户指定的第一零点
        first_null = abs(first_null_user);
    end
    
    % ========== 确定主瓣范围 ==========
    % 主瓣范围：[steer_angle - first_null, steer_angle + first_null]
    mainlobe_mask = (angles >= steer_angle - first_null - null_tolerance) & ...
                    (angles <= steer_angle + first_null + null_tolerance);
    
    % 旁瓣范围（主瓣范围外）
    sidelobe_mask = ~mainlobe_mask;
    
    % ========== 计算PSLR ==========
    if any(sidelobe_mask)
        % 找到旁瓣区域的最高值
        sidelobe_db = beam_pattern(sidelobe_mask);
        max_sidelobe_db = max(sidelobe_db);
        
        % PSLR = 主瓣峰值 - 最高旁瓣峰值
        PSLR = main_peak_db - max_sidelobe_db;
    else
        % 如果没有旁瓣区域，返回NaN
        warning('无法找到旁瓣区域，可能主瓣范围设置过大');
        PSLR = NaN;
    end
    
end

