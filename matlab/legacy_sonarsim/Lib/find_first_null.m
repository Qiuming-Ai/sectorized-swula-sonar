function [left_null, right_null] = find_first_null(beam_pattern, steer_angle, angle_range)
% FIND_FIRST_NULL 寻找波束图的第一零点
%
% 输入参数:
%   beam_pattern - 波束图数据 [N×1] (dB值)
%   steer_angle  - 扫描角度 (度)，通常为主瓣指向角度
%   angle_range  - 角度范围 [min_angle, max_angle] (度)
%                  或角度向量 [N×1] (度)
%
% 输出参数:
%   left_null    - steer角度左边的第一零点角度 (度)
%   right_null   - steer角度右边的第一零点角度 (度)
%
% 说明:
%   第一零点定义为离开主瓣后出现的第一个局部最小值（通常是最低点）
%   如果未找到零点，返回 NaN
%
% 示例:
%   % 使用角度范围
%   beam_db = ...;  % 波束图数据
%   [left, right] = find_first_null(beam_db, 0, [-85, 85]);
%
%   % 使用角度向量
%   angles = -85:0.1:85;
%   beam_db = ...;  % 波束图数据
%   [left, right] = find_first_null(beam_db, 0, angles);

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
    
    peak_val = beam_pattern(peak_idx);
    
    % ========== 找到第一零点（左右各一个） ==========
    % 策略：从峰值位置开始，找到离开主瓣后的第一个局部最小值
    
    left_null = NaN;
    right_null = NaN;
    
    % ========== 向左搜索第一零点 ==========
    if peak_idx > 2
        % 从峰值位置向左搜索第一个局部最小值
        % 跳过主瓣区域：找到值下降超过一定阈值后的区域
        
        % 从峰值向左，找到第一个低于峰值-10dB的点（大致离开主瓣）
        threshold_mainlobe = peak_val - 10;
        search_start_idx = peak_idx;
        for i = peak_idx:-1:1
            if beam_pattern(i) < threshold_mainlobe
                search_start_idx = i;
                break;
            end
        end
        
        % 从search_start_idx继续向左，搜索局部最小值
        % 局部最小值：相邻点都更大
        for i = search_start_idx:-1:2
            % 检查是否为局部最小值
            if i > 1 && i < N
                if beam_pattern(i) < beam_pattern(i-1) && beam_pattern(i) < beam_pattern(i+1)
                    % 找到局部最小值，检查是否足够低
                    if beam_pattern(i) < peak_val - 15 || beam_pattern(i) < -30
                        left_null = angles(i);
                        break;
                    end
                end
            elseif i == 1
                if N > 1 && beam_pattern(i) < beam_pattern(i+1)
                    if beam_pattern(i) < peak_val - 15 || beam_pattern(i) < -30
                        left_null = angles(i);
                        break;
                    end
                end
            end
        end
        
        % 如果还没找到，在左侧区域找最低点作为备选
        if isnan(left_null) && search_start_idx > 1
            search_end_idx = max(1, search_start_idx - 100);  % 限制搜索范围
            [min_val, min_local_idx] = min(beam_pattern(search_end_idx:search_start_idx));
            min_idx = search_end_idx + min_local_idx - 1;
            if min_val < peak_val - 15 && min_val < -25
                left_null = angles(min_idx);
            end
        end
    end
    
    % ========== 向右搜索第一零点 ==========
    if peak_idx < N-1
        % 从峰值位置向右搜索第一个局部最小值
        
        % 从峰值向右，找到第一个低于峰值-10dB的点（大致离开主瓣）
        threshold_mainlobe = peak_val - 10;
        search_start_idx = peak_idx;
        for i = peak_idx:N
            if beam_pattern(i) < threshold_mainlobe
                search_start_idx = i;
                break;
            end
        end
        
        % 从search_start_idx继续向右，搜索局部最小值
        for i = search_start_idx:N-1
            % 检查是否为局部最小值
            if i > 1 && i < N
                if beam_pattern(i) < beam_pattern(i-1) && beam_pattern(i) < beam_pattern(i+1)
                    % 找到局部最小值，检查是否足够低
                    if beam_pattern(i) < peak_val - 15 || beam_pattern(i) < -30
                        right_null = angles(i);
                        break;
                    end
                end
            elseif i == N
                if beam_pattern(i) < beam_pattern(i-1)
                    if beam_pattern(i) < peak_val - 15 || beam_pattern(i) < -30
                        right_null = angles(i);
                        break;
                    end
                end
            end
        end
        
        % 如果还没找到，在右侧区域找最低点作为备选
        if isnan(right_null) && search_start_idx < N
            search_end_idx = min(N, search_start_idx + 100);  % 限制搜索范围
            [min_val, min_local_idx] = min(beam_pattern(search_start_idx:search_end_idx));
            min_idx = search_start_idx + min_local_idx - 1;
            if min_val < peak_val - 15 && min_val < -25
                right_null = angles(min_idx);
            end
        end
    end
    
    % ========== 验证结果 ==========
    % 确保找到的零点确实在steer角度的左右两侧
    if ~isnan(left_null) && left_null > steer_angle
        left_null = NaN;
    end
    
    if ~isnan(right_null) && right_null < steer_angle
        right_null = NaN;
    end
    
end
