function plot_paperA_exp5_pulsecompression(varargin)
% PLOT_PAPERA_EXP5_PULSECOMPRESSION 读取并绘制 PaperA_exp5 脉冲压缩后的时域波形（科研绘图标准）
%
% 函数会自动读取 PaperA_exp5_0_pc.mat
% 包含三种编码方式的脉冲压缩结果：
%   - tw_db0_d: Baseline 100kHz LFM
%   - tw_db1_d: CDM 100kHz Zadoff-Chu编码
%   - tw_db2_d: FDM 20kHz LFM
%   - t_tw: 时间序列
%
% 可选参数 (通过 varargin 或结构体):
%   'dyn_range'    - 动态范围 (dB), 默认 60
%   'figure_handle' - 图形句柄 (可选，用于在指定图形中绘制)
%   'legend_labels' - 图例标签，默认 {'Baseline 100kHz LFM', 'CDM 100kHz Zadoff-Chu', 'FDM 20kHz LFM'}
%   'line_width'   - 线宽, 默认 2.0 (科研绘图推荐)
%   'font_size'     - 字体大小, 默认 14
%   'save_fig'      - 是否保存图片, 默认 false
%   'save_path'     - 保存路径, 默认 'pulsecompression_paperA_exp5.png'
%   'fig_size'      - 图形尺寸 [width, height] (inches), 默认 [10, 7]
%   'plot_type'     - 绘图类型 ('overlay' 或 'subplot'), 默认 'overlay'
%
% 示例:
%   plot_paperA_exp5_pulsecompression();
%   plot_paperA_exp5_pulsecompression('dyn_range', 40, 'save_fig', true, 'plot_type', 'subplot');

    % ========== 参数解析 ==========
    p = inputParser;
    addParameter(p, 'dyn_range', 60, @isnumeric);
    addParameter(p, 'figure_handle', [], @(x) isempty(x) || isgraphics(x));
    addParameter(p, 'legend_labels', {'Baseline 100kHz LFM', 'CDM 100kHz Zadoff-Chu', 'FDM 20kHz LFM'}, ...
        @(x) iscell(x) && length(x) == 3);
    addParameter(p, 'line_width', 2.0, @isnumeric);
    addParameter(p, 'font_size', 14, @isnumeric);
    addParameter(p, 'save_fig', false, @islogical);
    addParameter(p, 'save_path', 'pulsecompression_paperA_exp5.png', @(x) ischar(x) || isstring(x));
    addParameter(p, 'fig_size', [10, 7], @(x) isnumeric(x) && length(x) == 2);
    addParameter(p, 'plot_type', 'overlay', @(x) ismember(x, {'overlay', 'subplot'}));
    parse(p, varargin{:});
    
    dyn_range = p.Results.dyn_range;
    fig_handle = p.Results.figure_handle;
    legend_labels = p.Results.legend_labels;
    line_width = p.Results.line_width;
    font_size = p.Results.font_size;
    save_fig = p.Results.save_fig;
    save_path = p.Results.save_path;
    fig_size = p.Results.fig_size;
    plot_type = p.Results.plot_type;
    
    % ========== 默认文件路径 ==========
    mat_file = legacy_figure_asset('PaperA_exp5_0_pc.mat');
    
    % ========== 读取数据 ==========
    if ~exist(mat_file, 'file')
        error('文件不存在: %s', mat_file);
    end
    
    data = load(mat_file);
    
    % 检查必要的变量是否存在
    required_vars = {'tw_db0_d', 'tw_db1_d', 'tw_db2_d', 't_tw'};
    missing_vars = {};
    for i = 1:length(required_vars)
        if ~isfield(data, required_vars{i})
            missing_vars{end+1} = required_vars{i};
        end
    end
    
    if ~isempty(missing_vars)
        error('MAT文件中缺少以下变量: %s', strjoin(missing_vars, ', '));
    end
    
    % 提取数据
    tw_db0_d = data.tw_db0_d;  % Baseline 100kHz LFM
    tw_db1_d = data.tw_db1_d;  % CDM 100kHz Zadoff-Chu
    tw_db2_d = data.tw_db2_d;  % FDM 20kHz LFM
    t_tw = data.t_tw;          % 时间序列
    
    % 确保数据是列向量
    if size(t_tw, 1) == 1
        t_tw = t_tw(:);
    end
    if size(tw_db0_d, 1) == 1, tw_db0_d = tw_db0_d(:); end
    if size(tw_db1_d, 1) == 1, tw_db1_d = tw_db1_d(:); end
    if size(tw_db2_d, 1) == 1, tw_db2_d = tw_db2_d(:); end
    
    % 检查数据维度
    N = length(t_tw);
    if length(tw_db0_d) ~= N || length(tw_db1_d) ~= N || length(tw_db2_d) ~= N
        error('数据维度不匹配: t_tw长度=%d, 各tw_db_d长度不一致', N);
    end
    
    % 限制动态范围
    tw_db0_plot = max(tw_db0_d, -dyn_range);
    tw_db1_plot = max(tw_db1_d, -dyn_range);
    tw_db2_plot = max(tw_db2_d, -dyn_range);
    
    % ========== 计算PSLR、HPBW和ISLR ==========
    % 存储所有脉冲压缩数据
    pc_data = {tw_db0_d, tw_db1_d, tw_db2_d};
    pc_names = legend_labels;
    
    % 计算每个编码方式的PSLR、HPBW和ISLR
    % 其中 HPBW 换算为距离（单位: m），假设声速 c = 1500 m/s，采用双程传播：R = c * Δt / 2
    fprintf('\n========== 脉冲压缩PSLR、HPBW(距离)和ISLR计算结果 ==========\n');
    fprintf('%-30s | %10s | %12s | %10s\n', 'Method', 'PSLR(dB)', 'HPBW(m)', 'ISLR(dB)');
    fprintf('--------------------------------------------------------------------------\n');
    
    PSLR_results = zeros(3, 1);
    HPBW_results = zeros(3, 1);
    ISLR_results = zeros(3, 1);
    
    for i = 1:3
        pc_db = pc_data{i};
        
        % 找到主瓣峰值位置（通常在时间0附近）
        % 在时间0附近搜索峰值（±10%的时间范围）
        search_range_ratio = 0.1;  % 搜索范围：±10%的时间范围
        time_range = max(t_tw) - min(t_tw);
        search_range = search_range_ratio * time_range;
        peak_search_mask = abs(t_tw) <= search_range;
        peak_search_indices = find(peak_search_mask);
        
        if isempty(peak_search_indices)
            % 如果搜索范围太小，在整个范围内搜索
            [~, peak_idx] = max(pc_db);
        else
            [~, local_max_idx] = max(pc_db(peak_search_indices));
            peak_idx = peak_search_indices(local_max_idx);
        end
        
        main_peak_db = pc_db(peak_idx);
        peak_time = t_tw(peak_idx);
        
        % ========== 找到第一零点位置（用于定义主瓣范围） ==========
        % 方法1：寻找第一个局部最小值（更稳健）
        % 向左搜索第一个局部最小值
        left_null_idx = NaN;
        for j = peak_idx:-1:3
            % 检查是否为局部最小值（比前后两个点都小）
            if pc_db(j) < pc_db(j-1) && pc_db(j) < pc_db(j+1)
                left_null_idx = j;
                break;
            end
        end
        
        % 向右搜索第一个局部最小值
        right_null_idx = NaN;
        for j = peak_idx:(length(pc_db)-2)
            % 检查是否为局部最小值
            if pc_db(j) < pc_db(j-1) && pc_db(j) < pc_db(j+1)
                right_null_idx = j;
                break;
            end
        end
        
        % 方法2：如果找不到局部最小值，使用阈值方法（相对于峰值-30dB）
        threshold_null = main_peak_db - 30;  % 第一零点阈值
        
        if isnan(left_null_idx)
            % 向左搜索第一个低于阈值的点
            left_idx = find(pc_db(1:peak_idx) <= threshold_null, 1, 'last');
            if ~isempty(left_idx) && left_idx > 1
                % 使用线性插值
                if abs(pc_db(left_idx-1) - pc_db(left_idx)) > eps
                    t_interp = (threshold_null - pc_db(left_idx)) / (pc_db(left_idx-1) - pc_db(left_idx));
                    left_null_time = t_tw(left_idx) + t_interp * (t_tw(left_idx-1) - t_tw(left_idx));
                else
                    left_null_time = t_tw(left_idx);
                end
            elseif ~isempty(left_idx)
                left_null_time = t_tw(left_idx);
            else
                % 使用默认估计：峰值附近下降到-20dB的点
                threshold_estimate = main_peak_db - 20;
                left_idx = find(pc_db(1:peak_idx) <= threshold_estimate, 1, 'last');
                if ~isempty(left_idx)
                    left_null_time = t_tw(left_idx);
                else
                    left_null_time = t_tw(1);
                end
            end
        else
            left_null_time = t_tw(left_null_idx);
        end
        
        if isnan(right_null_idx)
            % 向右搜索第一个低于阈值的点
            right_idx = find(pc_db(peak_idx:end) <= threshold_null, 1, 'first');
            if ~isempty(right_idx)
                right_idx = peak_idx + right_idx - 1;
                if right_idx < length(pc_db) && abs(pc_db(right_idx+1) - pc_db(right_idx)) > eps
                    % 使用线性插值
                    t_interp = (threshold_null - pc_db(right_idx)) / (pc_db(right_idx+1) - pc_db(right_idx));
                    right_null_time = t_tw(right_idx) + t_interp * (t_tw(right_idx+1) - t_tw(right_idx));
                else
                    right_null_time = t_tw(right_idx);
                end
            else
                % 使用默认估计：峰值附近下降到-20dB的点
                threshold_estimate = main_peak_db - 20;
                right_idx = find(pc_db(peak_idx:end) <= threshold_estimate, 1, 'first');
                if ~isempty(right_idx)
                    right_null_time = t_tw(peak_idx + right_idx - 1);
                else
                    right_null_time = t_tw(end);
                end
            end
        else
            right_null_time = t_tw(right_null_idx);
        end
        
        % 计算第一零点距离（用于定义主瓣范围）
        dist_left = abs(left_null_time - peak_time);
        dist_right = abs(right_null_time - peak_time);
        first_null = max(dist_left, dist_right);
        
        % 定义主瓣范围
        mainlobe_mask = (t_tw >= peak_time - first_null) & (t_tw <= peak_time + first_null);
        sidelobe_mask = ~mainlobe_mask;
        
        % ========== 计算PSLR ==========
        if any(sidelobe_mask)
            sidelobe_db = pc_db(sidelobe_mask);
            max_sidelobe_db = max(sidelobe_db);
            PSLR = main_peak_db - max_sidelobe_db;
        else
            PSLR = NaN;
        end
        
        % ========== 计算HPBW ==========
        % 在主瓣范围内搜索-3dB点
        mainlobe_indices = find(mainlobe_mask);
        mainlobe_time = t_tw(mainlobe_indices);
        mainlobe_pc = pc_db(mainlobe_indices);
        
        % 确保时间是排序的
        [mainlobe_time, sort_idx] = sort(mainlobe_time);
        mainlobe_pc = mainlobe_pc(sort_idx);
        
        % 找到峰值在主瓣范围内的索引
        [~, peak_idx_in_mainlobe] = max(mainlobe_pc);
        
        % 计算-3dB阈值
        threshold_3dB = main_peak_db - 3.0;
        
        % 找到左侧-3dB点
        left_3dB_time = NaN;
        if peak_idx_in_mainlobe > 1
            for j = peak_idx_in_mainlobe:-1:2
                if mainlobe_pc(j) >= threshold_3dB && mainlobe_pc(j-1) < threshold_3dB
                    if abs(mainlobe_pc(j-1) - mainlobe_pc(j)) > eps
                        t = (threshold_3dB - mainlobe_pc(j)) / (mainlobe_pc(j-1) - mainlobe_pc(j));
                        left_3dB_time = mainlobe_time(j) + t * (mainlobe_time(j-1) - mainlobe_time(j));
                        break;
                    end
                elseif mainlobe_pc(j-1) <= threshold_3dB && mainlobe_pc(j) > threshold_3dB
                    if abs(mainlobe_pc(j-1) - mainlobe_pc(j)) > eps
                        t = (threshold_3dB - mainlobe_pc(j)) / (mainlobe_pc(j-1) - mainlobe_pc(j));
                        left_3dB_time = mainlobe_time(j) + t * (mainlobe_time(j-1) - mainlobe_time(j));
                        break;
                    end
                end
            end
            if isnan(left_3dB_time) && mainlobe_pc(1) <= threshold_3dB
                left_3dB_time = mainlobe_time(1);
            end
        end
        
        % 找到右侧-3dB点
        right_3dB_time = NaN;
        if peak_idx_in_mainlobe < length(mainlobe_pc)
            for j = peak_idx_in_mainlobe:(length(mainlobe_pc)-1)
                if mainlobe_pc(j) >= threshold_3dB && mainlobe_pc(j+1) < threshold_3dB
                    if abs(mainlobe_pc(j+1) - mainlobe_pc(j)) > eps
                        t = (threshold_3dB - mainlobe_pc(j)) / (mainlobe_pc(j+1) - mainlobe_pc(j));
                        right_3dB_time = mainlobe_time(j) + t * (mainlobe_time(j+1) - mainlobe_time(j));
                        break;
                    end
                elseif mainlobe_pc(j+1) <= threshold_3dB && mainlobe_pc(j) > threshold_3dB
                    if abs(mainlobe_pc(j+1) - mainlobe_pc(j)) > eps
                        t = (threshold_3dB - mainlobe_pc(j)) / (mainlobe_pc(j+1) - mainlobe_pc(j));
                        right_3dB_time = mainlobe_time(j) + t * (mainlobe_time(j+1) - mainlobe_time(j));
                        break;
                    end
                end
            end
            if isnan(right_3dB_time) && mainlobe_pc(end) <= threshold_3dB
                right_3dB_time = mainlobe_time(end);
            end
        end
        
        % 计算HPBW（时间宽度换算为距离，声速 c = 1500 m/s，双程：R = c * Δt / 2）
        if ~isnan(left_3dB_time) && ~isnan(right_3dB_time)
            delta_t = (right_3dB_time - left_3dB_time);  % 单位: s
            c_sound = 1500;                              % 声速 (m/s)
            HPBW = c_sound * delta_t / 2;                % 转换为距离 (m)，双程传播
        else
            HPBW = NaN;
        end
        
        % ========== 计算ISLR ==========
        % 使用功率积分：power = 10^(dB/10)，然后分别在主瓣区和旁瓣区求和
        if any(mainlobe_mask) && any(sidelobe_mask)
            pc_power = 10.^(pc_db/10);                % 相对功率
            main_power = sum(pc_power(mainlobe_mask));
            side_power = sum(pc_power(sidelobe_mask));
            if main_power > 0 && side_power > 0
                ISLR = 10*log10(side_power / main_power);
            else
                ISLR = NaN;
            end
        else
            ISLR = NaN;
        end
        
        PSLR_results(i) = PSLR;
        HPBW_results(i) = HPBW;
        ISLR_results(i) = ISLR;
        
        fprintf('%-30s | %10.2f | %12.3f | %10.2f\n', ...
                pc_names{i}, PSLR, HPBW, ISLR);
    end
    fprintf('====================================================\n\n');
    
    % ========== 科研绘图颜色方案 ==========
    % 使用科研论文常用的颜色方案（3种颜色）
    colors = [0.2, 0.4, 0.8;      % 蓝色 - Baseline
              0.8, 0.2, 0.2;      % 红色 - CDM
              0.2, 0.6, 0.3];     % 绿色 - FDM
    
    % 线型：实线、虚线、点划线
    line_styles = {'-', '--', '-.'};
    
    % ========== 绘图 ==========
    if isempty(fig_handle)
        fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, fig_size(1), fig_size(2)]);
    else
        fig = figure(fig_handle);
    end
    
    if strcmp(plot_type, 'subplot')
        % 子图模式：每个编码方式单独一个子图
        for i = 1:3
            ax = subplot(3, 1, i);
            
            switch i
                case 1
                    plot(ax, t_tw*1e6, tw_db0_plot, 'LineWidth', line_width, ...
                        'Color', colors(1,:), 'LineStyle', line_styles{1});
                    title_str = legend_labels{1};
                case 2
                    plot(ax, t_tw*1e6, tw_db1_plot, 'LineWidth', line_width, ...
                        'Color', colors(2,:), 'LineStyle', line_styles{2});
                    title_str = legend_labels{2};
                case 3
                    plot(ax, t_tw*1e6, tw_db2_plot, 'LineWidth', line_width, ...
                        'Color', colors(3,:), 'LineStyle', line_styles{3});
                    title_str = legend_labels{3};
            end
            
            % 坐标轴标签和标题
            if i == 3
                xlabel('Time (\mus)', 'FontSize', font_size, 'FontWeight', 'bold');
            end
            ylabel('Amplitude (dB)', 'FontSize', font_size, 'FontWeight', 'bold');
            title(title_str, 'FontSize', font_size, 'FontWeight', 'bold');
            
            % 网格设置
            grid(ax, 'on');
            grid(ax, 'minor');
            ax.GridAlpha = 0.3;
            ax.MinorGridAlpha = 0.15;
            ax.GridLineStyle = '-';
            ax.MinorGridLineStyle = ':';
            
            % 坐标轴范围
            ylim(ax, [-dyn_range, 5]);
            xlim(ax, [min(t_tw)*1e6, max(t_tw)*1e6]);
            
            % 坐标轴字体大小
            ax.FontSize = font_size-1;
            ax.TickLength = [0.02, 0.02];
            
            % 边框设置
            box(ax, 'on');
            
            % 添加零分贝参考线
            if exist('yline', 'builtin') || exist('yline', 'file')
                yline(ax, 0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8, 'Alpha', 0.5, ...
                    'HandleVisibility', 'off');
            else
                plot(ax, [min(t_tw)*1e6, max(t_tw)*1e6], [0, 0], '--', ...
                    'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
        end
        
        % 整体标题
        sgtitle('Pulse Compression Results Comparison', 'FontSize', font_size+2, 'FontWeight', 'bold');
        
    else
        % 叠加模式：所有三条线在同一图中
        ax = axes('Parent', fig);
        
        h1 = plot(ax, t_tw*1e6, tw_db0_plot, 'LineWidth', line_width, ...
            'Color', colors(1,:), 'LineStyle', line_styles{1}); hold on;
        h2 = plot(ax, t_tw*1e6, tw_db1_plot, 'LineWidth', line_width, ...
            'Color', colors(2,:), 'LineStyle', line_styles{2});
        h3 = plot(ax, t_tw*1e6, tw_db2_plot, 'LineWidth', line_width, ...
            'Color', colors(3,:), 'LineStyle', line_styles{3});
        
        % 图例设置（在绘制参考线之前创建）
        legend(ax, [h1, h2, h3], legend_labels, 'Location', 'best', ...
            'FontSize', font_size-1, 'Box', 'on', 'EdgeColor', [0.3, 0.3, 0.3]);
        
        % 添加零分贝参考线（兼容旧版本MATLAB，不显示在图例中）
        if exist('yline', 'builtin') || exist('yline', 'file')
            yline(ax, 0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8, 'Alpha', 0.5, ...
                'HandleVisibility', 'off');
        else
            plot(ax, [min(t_tw)*1e6, max(t_tw)*1e6], [0, 0], '--', ...
                'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        % 坐标轴标签和标题
        xlabel('Time (\mus)', 'FontSize', font_size, 'FontWeight', 'bold');
        ylabel('Amplitude (dB)', 'FontSize', font_size, 'FontWeight', 'bold');
        title('Pulse Compression Results Comparison', 'FontSize', font_size+1, 'FontWeight', 'bold');
        
        % 网格设置（科研绘图常用细网格）
        grid(ax, 'on');
        grid(ax, 'minor');
        ax.GridAlpha = 0.3;
        ax.MinorGridAlpha = 0.15;
        ax.GridLineStyle = '-';
        ax.MinorGridLineStyle = ':';
        
        % 坐标轴范围
        ylim(ax, [-dyn_range, 5]);
        xlim(ax, [min(t_tw)*1e6, max(t_tw)*1e6]);
        
        % 坐标轴字体大小
        ax.FontSize = font_size-1;
        ax.TickLength = [0.02, 0.02];
        
        % 边框设置
        box(ax, 'on');
        
        hold off;
    end
    
    % 设置图形背景
    set(fig, 'Color', 'w');
    
    % ========== 保存图片 ==========
    if save_fig
        % 设置保存选项（高质量）
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, save_path, '-dpng', '-r300');  % 300 DPI 高质量
        fprintf('图片已保存至: %s\n', save_path);
    end
    
end

