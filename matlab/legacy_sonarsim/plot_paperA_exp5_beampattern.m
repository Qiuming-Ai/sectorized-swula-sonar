function plot_paperA_exp5_beampattern(varargin)
% PLOT_PAPERA_EXP5_BEAMPATTERN 读取并绘制 PaperA_exp5 的波束图（科研绘图标准）
%
% 函数会自动读取 PaperA_exp5_0.mat 和 PaperA_exp5_1.mat
%
% 可选参数 (通过 varargin 或结构体):
%   'dyn_range'    - 动态范围 (dB), 默认 60
%   'figure_handle' - 图形句柄 (可选，用于在指定图形中绘制)
%   'legend_labels' - 图例标签，默认 {'Baseline', 'CDM', 'FDM', 'Random', 'Nested', 'Coprime'}
%   'line_width'   - 线宽, 默认 2.0 (科研绘图推荐)
%   'plot_type'    - 绘图类型 ('rect' 或 'polar'), 默认 'rect'
%   'font_size'     - 字体大小, 默认 14
%   'save_fig'      - 是否保存图片, 默认 false
%   'save_path'     - 保存路径, 默认 'beampattern_comparison.png'
%   'fig_size'      - 图形尺寸 [width, height] (inches), 默认 [10, 7]
%
% 示例:
%   plot_paperA_exp5_beampattern();
%   plot_paperA_exp5_beampattern('PaperA_exp5_0.mat', 'dyn_range', 40, 'save_fig', true);

    % ========== 参数解析 ==========
    p = inputParser;
    addParameter(p, 'dyn_range', 60, @isnumeric);
    addParameter(p, 'figure_handle', [], @(x) isempty(x) || isgraphics(x));
    addParameter(p, 'legend_labels', {'Baseline', 'CDM', 'FDM', 'Random', 'Nested', 'Coprime'}, ...
        @(x) iscell(x) && length(x) == 6);
    addParameter(p, 'line_width', 2.0, @isnumeric);
    addParameter(p, 'plot_type', 'rect', @(x) ismember(x, {'polar', 'rect'}));
    addParameter(p, 'font_size', 14, @isnumeric);
    addParameter(p, 'save_fig', false, @islogical);
    addParameter(p, 'save_path', 'beampattern_comparison.png', @(x) ischar(x) || isstring(x));
    addParameter(p, 'fig_size', [10, 7], @(x) isnumeric(x) && length(x) == 2);
    parse(p, varargin{:});
    
    dyn_range = p.Results.dyn_range;
    fig_handle = p.Results.figure_handle;
    legend_labels = p.Results.legend_labels;
    line_width = p.Results.line_width;
    plot_type = p.Results.plot_type;
    font_size = p.Results.font_size;
    save_fig = p.Results.save_fig;
    save_path = p.Results.save_path;
    fig_size = p.Results.fig_size;
    
    % ========== 默认文件路径 ==========
    mat_file_0 = legacy_figure_asset('PaperA_exp5_0.mat');
    mat_file_1 = legacy_figure_asset('PaperA_exp5_1.mat');
    
    % ========== 读取第一个文件 (PaperA_exp5_0.mat) ==========
    if ~exist(mat_file_0, 'file')
        error('文件不存在: %s', mat_file_0);
    end
    
    data0 = load(mat_file_0);
    
    % 检查必要的变量是否存在
    required_vars_0 = {'tw_db0', 'tw_db1', 'tw_db2', 'theta_deg'};
    missing_vars = {};
    for i = 1:length(required_vars_0)
        if ~isfield(data0, required_vars_0{i})
            missing_vars{end+1} = required_vars_0{i};
        end
    end
    
    if ~isempty(missing_vars)
        error('MAT文件 %s 中缺少以下变量: %s', mat_file_0, strjoin(missing_vars, ', '));
    end
    
    % 提取第一个文件的数据
    tw_db0 = data0.tw_db0;
    tw_db1 = data0.tw_db1+4.8128;
    % tw_db1 = circshift(tw_db1,1);
    tw_db2 = data0.tw_db2+5.584;
    % tw_db2 = circshift(tw_db2,1);
    theta_deg = data0.theta_deg;
    
    % ========== 读取第二个文件 (PaperA_exp5_1.mat) ==========
    if ~exist(mat_file_1, 'file')
        error('文件不存在: %s', mat_file_1);
    end
    
    data1 = load(mat_file_1);
    
    % 检查必要的变量是否存在
    required_vars_1 = {'tw_db3', 'tw_db4', 'tw_db5'};
    missing_vars = {};
    for i = 1:length(required_vars_1)
        if ~isfield(data1, required_vars_1{i})
            missing_vars{end+1} = required_vars_1{i};
        end
    end
    
    if ~isempty(missing_vars)
        error('MAT文件 %s 中缺少以下变量: %s', mat_file_1, strjoin(missing_vars, ', '));
    end
    
    % 提取第二个文件的数据
    tw_db3 = data1.tw_db3;  % Random
    tw_db4 = data1.tw_db4;  % Nested
    tw_db5 = data1.tw_db5;  % Coprime
    
    % 确保数据是列向量
    if size(theta_deg, 1) == 1
        theta_deg = theta_deg(:);
    end
    if size(tw_db0, 1) == 1, tw_db0 = tw_db0(:); end
    if size(tw_db1, 1) == 1, tw_db1 = tw_db1(:); end
    if size(tw_db2, 1) == 1, tw_db2 = tw_db2(:); end
    if size(tw_db3, 1) == 1, tw_db3 = tw_db3(:); end
    if size(tw_db4, 1) == 1, tw_db4 = tw_db4(:); end
    if size(tw_db5, 1) == 1, tw_db5 = tw_db5(:); end
    
    % 检查数据维度
    N = length(theta_deg);
    if length(tw_db0) ~= N || length(tw_db1) ~= N || length(tw_db2) ~= N || ...
       length(tw_db3) ~= N || length(tw_db4) ~= N || length(tw_db5) ~= N
        error('数据维度不匹配: theta_deg长度=%d, 各tw_db长度不一致', N);
    end
    
    % ========== 计算PSLR和ISLR ==========
    % 第一零点位置（度）
    first_null = struct();
    first_null.Baseline = 1.8;
    first_null.CDM = 1.9;
    first_null.FDM = 1.9;
    first_null.Random = 1.1;
    first_null.Nested = 1.3;
    first_null.Coprime = 0.9;
    
    % 存储所有波束数据
    beam_data = {tw_db0, tw_db1, tw_db2, tw_db3, tw_db4, tw_db5};
    beam_names = {'Baseline', 'CDM', 'FDM', 'Random', 'Nested', 'Coprime'};
    % beam_data = {tw_db0, tw_db1};
    % beam_names = {'Baseline', 'CDM'};
    first_nulls = [first_null.Baseline, first_null.CDM, first_null.FDM, ...
                   first_null.Random, first_null.Nested, first_null.Coprime];
    
    % 计算每个波束的PSLR、ISLR和HPBW
    fprintf('\n========== PSLR、ISLR和HPBW计算结果 ==========\n');
    fprintf('%-12s | %8s | %8s | %8s | %12s\n', 'Method', 'PSLR(dB)', 'ISLR(dB)', 'HPBW(°)', 'First Null(°)');
    fprintf('--------------------------------------------------------\n');
    
    PSLR_results = zeros(6, 1);
    ISLR_results = zeros(6, 1);
    HPBW_results = zeros(6, 1);
    
    for i = 1:6
        beam_db = beam_data{i};
        first_null_deg = first_nulls(i);
        
        % 找到主瓣峰值位置（通常在0°附近）
        [~, main_peak_idx] = max(beam_db);
        main_peak_db = beam_db(main_peak_idx);
        
        % 确定主瓣范围（-first_null 到 +first_null）
        mainlobe_mask = (theta_deg >= -first_null_deg) & (theta_deg <= first_null_deg);
        
        % 旁瓣范围（主瓣范围外）
        sidelobe_mask = ~mainlobe_mask;
        
        % 计算PSLR：主瓣峰值与最高旁瓣峰值的比值（dB）
        if any(sidelobe_mask)
            sidelobe_db = beam_db(sidelobe_mask);
            max_sidelobe_db = max(sidelobe_db);
            PSLR = main_peak_db - max_sidelobe_db;
        else
            PSLR = NaN;
        end
        
        % 计算ISLR：主瓣能量与旁瓣能量的比值（dB）
        % 将dB转换为线性值计算能量
        beam_linear = 10.^(beam_db / 20);  % 转换为幅度
        beam_power = beam_linear.^2;       % 转换为功率
        
        % 计算主瓣能量和旁瓣能量
        mainlobe_energy = sum(beam_power(mainlobe_mask));
        sidelobe_energy = sum(beam_power(sidelobe_mask));
        
        if sidelobe_energy > 0
            ISLR = 10*log10(mainlobe_energy / sidelobe_energy);
        else
            ISLR = Inf;
        end
        
        % 计算HPBW：半功率波束宽度（-3dB点之间的角度宽度）
        % 从第一零点位置向内搜索-3dB点
        threshold_3dB = main_peak_db - 3.0;  % -3dB阈值
        
        % 在主瓣范围内搜索（-first_null 到 +first_null）
        mainlobe_indices = find(mainlobe_mask);
        mainlobe_theta = theta_deg(mainlobe_indices);
        mainlobe_beam = beam_db(mainlobe_indices);
        
        % 确保角度是排序的（从小到大）
        [mainlobe_theta, sort_idx] = sort(mainlobe_theta);
        mainlobe_beam = mainlobe_beam(sort_idx);
        
        % 找到峰值在主瓣范围内的索引
        [~, peak_idx_in_mainlobe] = max(mainlobe_beam);
        
        % 找到左侧-3dB点（从峰值向左搜索，即从峰值索引向前搜索）
        left_3dB_angle = NaN;
        if peak_idx_in_mainlobe > 1
            % 从峰值位置向左（向前）搜索
            for j = peak_idx_in_mainlobe:-1:2
                if mainlobe_beam(j) >= threshold_3dB && mainlobe_beam(j-1) < threshold_3dB
                    % 跨越-3dB点，使用线性插值
                    t = (threshold_3dB - mainlobe_beam(j)) / (mainlobe_beam(j-1) - mainlobe_beam(j) + eps);
                    left_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j-1) - mainlobe_theta(j));
                    break;
                elseif mainlobe_beam(j-1) <= threshold_3dB && mainlobe_beam(j) > threshold_3dB
                    % 另一种跨越情况
                    t = (threshold_3dB - mainlobe_beam(j)) / (mainlobe_beam(j-1) - mainlobe_beam(j) + eps);
                    left_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j-1) - mainlobe_theta(j));
                    break;
                end
            end
            % 如果没有找到跨越点，检查是否第一个点就低于阈值
            if isnan(left_3dB_angle) && mainlobe_beam(1) <= threshold_3dB
                left_3dB_angle = mainlobe_theta(1);
            end
        end
        
        % 找到右侧-3dB点（从峰值向右搜索，即从峰值索引向后搜索）
        right_3dB_angle = NaN;
        if peak_idx_in_mainlobe < length(mainlobe_beam)
            % 从峰值位置向右（向后）搜索
            for j = peak_idx_in_mainlobe:(length(mainlobe_beam)-1)
                if mainlobe_beam(j) >= threshold_3dB && mainlobe_beam(j+1) < threshold_3dB
                    % 跨越-3dB点，使用线性插值
                    t = (threshold_3dB - mainlobe_beam(j)) / (mainlobe_beam(j+1) - mainlobe_beam(j) + eps);
                    right_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j+1) - mainlobe_theta(j));
                    break;
                elseif mainlobe_beam(j+1) <= threshold_3dB && mainlobe_beam(j) > threshold_3dB
                    % 另一种跨越情况
                    t = (threshold_3dB - mainlobe_beam(j)) / (mainlobe_beam(j+1) - mainlobe_beam(j) + eps);
                    right_3dB_angle = mainlobe_theta(j) + t * (mainlobe_theta(j+1) - mainlobe_theta(j));
                    break;
                end
            end
            % 如果没有找到跨越点，检查是否最后一个点就低于阈值
            if isnan(right_3dB_angle) && mainlobe_beam(end) <= threshold_3dB
                right_3dB_angle = mainlobe_theta(end);
            end
        end
        
        % 计算HPBW
        if ~isnan(left_3dB_angle) && ~isnan(right_3dB_angle)
            HPBW = right_3dB_angle - left_3dB_angle;
        else
            HPBW = NaN;
        end
        
        PSLR_results(i) = PSLR;
        ISLR_results(i) = ISLR;
        HPBW_results(i) = HPBW;
        
        fprintf('%-12s | %8.2f | %8.2f | %8.2f | %12.1f\n', ...
                beam_names{i}, PSLR, ISLR, HPBW, first_null_deg);
    end
    fprintf('====================================================\n\n');
    
    % 限制动态范围
    tw_db0_plot = max(tw_db0, -dyn_range);
    tw_db1_plot = max(tw_db1, -dyn_range);
    tw_db2_plot = max(tw_db2, -dyn_range);
    tw_db3_plot = max(tw_db3, -dyn_range);
    tw_db4_plot = max(tw_db4, -dyn_range);
    tw_db5_plot = max(tw_db5, -dyn_range);
    
    % ========== 科研绘图颜色方案 ==========
    % 使用科研论文常用的颜色方案（6种颜色）
    colors = [0.2, 0.4, 0.8;      % 蓝色 - Baseline
              0.8, 0.2, 0.2;      % 红色 - CDM
              0.2, 0.6, 0.3;      % 绿色 - FDM
              0.9, 0.6, 0.1;      % 橙色 - Random
              0.6, 0.2, 0.8;      % 紫色 - Nested
              0.0, 0.7, 0.7];     % 青色 - Coprime
    
    % 线型：实线、虚线、点划线、点线等
    line_styles = {'-', '--', '-.', ':', '-', '--'};
    
    % ========== 绘图 ==========
    if isempty(fig_handle)
        fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, fig_size(1), fig_size(2)]);
    else
        fig = figure(fig_handle);
    end
    
    ax = axes('Parent', fig);
    
    if strcmp(plot_type, 'polar')
        % 极坐标图
        % 将dB转换为归一化幅度 (0-1)
        tw_norm0 = 10.^(tw_db0_plot / 20);
        tw_norm1 = 10.^(tw_db1_plot / 20);
        tw_norm2 = 10.^(tw_db2_plot / 20);
        tw_norm3 = 10.^(tw_db3_plot / 20);
        tw_norm4 = 10.^(tw_db4_plot / 20);
        tw_norm5 = 10.^(tw_db5_plot / 20);
        
        % 归一化到最大值
        max_val = max([max(tw_norm0), max(tw_norm1), max(tw_norm2), ...
                      max(tw_norm3), max(tw_norm4), max(tw_norm5)]);
        tw_norm0 = tw_norm0 / max_val;
        tw_norm1 = tw_norm1 / max_val;
        tw_norm2 = tw_norm2 / max_val;
        tw_norm3 = tw_norm3 / max_val;
        tw_norm4 = tw_norm4 / max_val;
        tw_norm5 = tw_norm5 / max_val;
        
        % 兼容性检查：使用 polarplot (R2016a+) 或 polar (旧版本)
        if exist('polarplot', 'builtin') || exist('polarplot', 'file')
            polarplot(ax, deg2rad(theta_deg), tw_norm0, 'LineWidth', line_width, ...
                'Color', colors(1,:), 'LineStyle', line_styles{1}); hold on;
            polarplot(ax, deg2rad(theta_deg), tw_norm1, 'LineWidth', line_width, ...
                'Color', colors(2,:), 'LineStyle', line_styles{2});
            polarplot(ax, deg2rad(theta_deg), tw_norm2, 'LineWidth', line_width, ...
                'Color', colors(3,:), 'LineStyle', line_styles{3});
            polarplot(ax, deg2rad(theta_deg), tw_norm3, 'LineWidth', line_width, ...
                'Color', colors(4,:), 'LineStyle', line_styles{4});
            polarplot(ax, deg2rad(theta_deg), tw_norm4, 'LineWidth', line_width, ...
                'Color', colors(5,:), 'LineStyle', line_styles{5});
            polarplot(ax, deg2rad(theta_deg), tw_norm5, 'LineWidth', line_width, ...
                'Color', colors(6,:), 'LineStyle', line_styles{6});
            rlim([0, 1]);
        else
            % 旧版本使用 polar 函数
            polar(deg2rad(theta_deg), tw_norm0, 'LineWidth', line_width); hold on;
            polar(deg2rad(theta_deg), tw_norm1, 'LineWidth', line_width);
            polar(deg2rad(theta_deg), tw_norm2, 'LineWidth', line_width);
            polar(deg2rad(theta_deg), tw_norm3, 'LineWidth', line_width);
            polar(deg2rad(theta_deg), tw_norm4, 'LineWidth', line_width);
            polar(deg2rad(theta_deg), tw_norm5, 'LineWidth', line_width);
        end
        title('Beam Pattern (Polar)', 'FontSize', font_size, 'FontWeight', 'bold');
        legend(legend_labels, 'Location', 'best', 'FontSize', font_size-1, 'Box', 'on');
        grid on;
        hold off;
        
    else
        % 直角坐标图 - 科研绘图标准格式（6条线）
        h1 = plot(ax, theta_deg, tw_db0_plot, 'LineWidth', line_width, ...
            'Color', colors(1,:), 'LineStyle', line_styles{1}); hold on;
        h2 = plot(ax, theta_deg, tw_db1_plot, 'LineWidth', line_width, ...
            'Color', colors(2,:), 'LineStyle', line_styles{2});
        h3 = plot(ax, theta_deg, tw_db2_plot, 'LineWidth', line_width, ...
            'Color', colors(3,:), 'LineStyle', line_styles{3});
        h4 = plot(ax, theta_deg, tw_db3_plot, 'LineWidth', line_width, ...
            'Color', colors(4,:), 'LineStyle', line_styles{4});
        h5 = plot(ax, theta_deg, tw_db4_plot, 'LineWidth', line_width, ...
            'Color', colors(5,:), 'LineStyle', line_styles{5});
        h6 = plot(ax, theta_deg, tw_db5_plot, 'LineWidth', line_width, ...
            'Color', colors(6,:), 'LineStyle', line_styles{6});
        
        % 图例设置（在绘制参考线之前创建，只包含6条数据线）
        legend(ax, [h1, h2, h3, h4, h5, h6], legend_labels, 'Location', 'best', ...
            'FontSize', font_size-1, 'Box', 'on', 'EdgeColor', [0.3, 0.3, 0.3]);
        
        % 添加零分贝参考线（兼容旧版本MATLAB，不显示在图例中）
        if exist('yline', 'builtin') || exist('yline', 'file')
            yline(ax, 0, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8, 'Alpha', 0.5, ...
                'HandleVisibility', 'off');
        else
            % 旧版本使用 plot 绘制参考线，设置 HandleVisibility 为 off 避免出现在图例中
            plot(ax, [min(theta_deg), max(theta_deg)], [0, 0], '--', ...
                'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        
        % 坐标轴标签和标题
        xlabel('Angle (degrees)', 'FontSize', font_size, 'FontWeight', 'bold');
        ylabel('Beam Response (dB)', 'FontSize', font_size, 'FontWeight', 'bold');
        title('Beam Pattern Comparison', 'FontSize', font_size+1, 'FontWeight', 'bold');
        
        % 网格设置（科研绘图常用细网格）
        grid(ax, 'on');
        grid(ax, 'minor');
        ax.GridAlpha = 0.3;
        ax.MinorGridAlpha = 0.15;
        ax.GridLineStyle = '-';
        ax.MinorGridLineStyle = ':';
        
        % 坐标轴范围
        ylim(ax, [-dyn_range, 5]);
        xlim(ax, [min(theta_deg), max(theta_deg)]);
        
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

