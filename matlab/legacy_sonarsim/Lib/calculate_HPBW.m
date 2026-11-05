function HPBW = calculate_HPBW(beam_pattern, steer_angle, angle_range, varargin)
% CALCULATE_HPBW 计算波束图的半功率波束宽度（Half Power Beam Width）
%
% 输入参数:
%   beam_pattern - 波束图数据 [N×1] (dB值)
%   steer_angle  - 扫描角度 (度)，通常为主瓣指向角度
%   angle_range  - 角度范围 [min_angle, max_angle] (度)
%                  或角度向量 [N×1] (度)
%
% 可选参数:
%   'first_null' - 第一零点位置 (度)，如果提供则直接使用来定义主瓣范围
%                  如果不提供，将自动调用 find_first_null 函数计算
%   'null_tolerance' - 第一零点容差 (度)，用于定义主瓣范围，默认使用第一零点
%   'threshold_db' - 功率下降阈值 (dB)，默认 -3.0 (半功率点)
%
% 输出参数:
%   HPBW         - 半功率波束宽度 (度)，即-3dB点之间的角度宽度
%
% 说明:
%   HPBW (Half Power Beam Width) 定义为：
%   在主瓣内，从峰值下降到-3dB（半功率点）的两个角度之间的宽度
%   使用线性插值来精确计算-3dB点的位置
%
% 示例:
%   % 使用角度范围，自动计算第一零点
%   beam_db = ...;  % 波束图数据
%   HPBW = calculate_HPBW(beam_db, 0, [-85, 85]);
%
%   % 使用角度向量，自动计算第一零点
%   angles = -85:0.1:85;
%   beam_db = ...;  % 波束图数据
%   HPBW = calculate_HPBW(beam_db, 0, angles);
%
%   % 指定第一零点位置
%   HPBW = calculate_HPBW(beam_db, 0, angles, 'first_null', 1.8);

    % ========== 参数解析 ==========
    p = inputParser;
    addParameter(p, 'first_null', [], @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'null_tolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'threshold_db', -3.0, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});
    
    first_null_user = p.Results.first_null;
    null_tolerance = p.Results.null_tolerance;
    threshold_db = p.Results.threshold_db;
    
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
    peak_angle = angles(peak_idx);
    
    % ========== 确定第一零点位置（用于定义主瓣范围） ==========
    if isempty(first_null_user)
        % 自动计算第一零点
        try
            [left_null, right_null] = find_first_null(beam_pattern, steer_angle, angle_range);
            
            % 计算从steer角度到左右零点的距离，取较大者作为第一零点
            if ~isnan(left_null) && ~isnan(right_null)
                dist_left = abs(left_null - steer_angle);
                dist_right = abs(right_null - steer_angle);
                first_null = max(dist_left, dist_right);
            elseif ~isnan(left_null)
                first_null = abs(left_null - steer_angle);
            elseif ~isnan(right_null)
                first_null = abs(right_null - steer_angle);
            else
                % 如果找不到第一零点，使用默认值
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
    
    % 在主瓣范围内搜索-3dB点
    mainlobe_indices = find(mainlobe_mask);
    mainlobe_theta = angles(mainlobe_indices);
    mainlobe_beam = beam_pattern(mainlobe_indices);
    
    % 确保角度是排序的（从小到大）
    [mainlobe_theta, sort_idx] = sort(mainlobe_theta);
    mainlobe_beam = mainlobe_beam(sort_idx);
    
    % 找到峰值在主瓣范围内的索引
    [~, peak_idx_in_mainlobe] = max(mainlobe_beam);
    
    % ========== 计算阈值 ==========
    threshold_value = main_peak_db + threshold_db;  % -3dB阈值
    
    % ========== 找到左侧-3dB点 ==========
    left_3dB_angle = NaN;
    if peak_idx_in_mainlobe > 1
        % 从峰值位置向左搜索（从峰值索引向前搜索）
        for j = peak_idx_in_mainlobe:-1:2
            % 检查是否跨越阈值
            if mainlobe_beam(j) >= threshold_value && mainlobe_beam(j-1) < threshold_value
                % 跨越阈值点，使用线性插值
                if abs(mainlobe_beam(j-1) - mainlobe_beam(j)) > eps
                    t = (threshold_value - mainlobe_beam(j)) / (mainlobe_beam(j-1) - mainlobe_beam(j));
                    left_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j-1) - mainlobe_theta(j));
                    break;
                else
                    % 如果两个点值相同，取中点
                    left_3dB_angle = (mainlobe_theta(j) + mainlobe_theta(j-1)) / 2;
                    break;
                end
            elseif mainlobe_beam(j-1) <= threshold_value && mainlobe_beam(j) > threshold_value
                % 另一种跨越情况
                if abs(mainlobe_beam(j-1) - mainlobe_beam(j)) > eps
                    t = (threshold_value - mainlobe_beam(j)) / (mainlobe_beam(j-1) - mainlobe_beam(j));
                    left_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j-1) - mainlobe_theta(j));
                    break;
                else
                    left_3dB_angle = (mainlobe_theta(j) + mainlobe_theta(j-1)) / 2;
                    break;
                end
            end
        end
        
        % 如果没有找到跨越点，检查是否第一个点就低于阈值
        if isnan(left_3dB_angle) && mainlobe_beam(1) <= threshold_value
            left_3dB_angle = mainlobe_theta(1);
        end
    end
    
    % ========== 找到右侧-3dB点 ==========
    right_3dB_angle = NaN;
    if peak_idx_in_mainlobe < length(mainlobe_beam)
        % 从峰值位置向右搜索（从峰值索引向后搜索）
        for j = peak_idx_in_mainlobe:(length(mainlobe_beam)-1)
            % 检查是否跨越阈值
            if mainlobe_beam(j) >= threshold_value && mainlobe_beam(j+1) < threshold_value
                % 跨越阈值点，使用线性插值
                if abs(mainlobe_beam(j+1) - mainlobe_beam(j)) > eps
                    t = (threshold_value - mainlobe_beam(j)) / (mainlobe_beam(j+1) - mainlobe_beam(j));
                    right_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j+1) - mainlobe_theta(j));
                    break;
                else
                    % 如果两个点值相同，取中点
                    right_3dB_angle = (mainlobe_theta(j) + mainlobe_theta(j+1)) / 2;
                    break;
                end
            elseif mainlobe_beam(j+1) <= threshold_value && mainlobe_beam(j) > threshold_value
                % 另一种跨越情况
                if abs(mainlobe_beam(j+1) - mainlobe_beam(j)) > eps
                    t = (threshold_value - mainlobe_beam(j)) / (mainlobe_beam(j+1) - mainlobe_beam(j));
                    right_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j+1) - mainlobe_theta(j));
                    break;
                else
                    right_3dB_angle = (mainlobe_theta(j) + mainlobe_theta(j+1)) / 2;
                    break;
                end
            end
        end
        
        % 如果没有找到跨越点，检查是否最后一个点就低于阈值
        if isnan(right_3dB_angle) && mainlobe_beam(end) <= threshold_value
            right_3dB_angle = mainlobe_theta(end);
        end
    end
    
    % ========== 计算HPBW ==========
    if ~isnan(left_3dB_angle) && ~isnan(right_3dB_angle)
        HPBW = right_3dB_angle - left_3dB_angle;
    else
        % 如果找不到-3dB点，返回NaN
        warning('无法找到完整的-3dB点，HPBW计算结果为NaN');
        HPBW = NaN;
    end
    
end










