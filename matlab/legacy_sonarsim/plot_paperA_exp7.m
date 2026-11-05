function plot_paperA_exp7(varargin)
% PLOT_PAPERA_EXP7 读取并绘制 PaperA_exp7_0.mat 的PSLR和HPBW随扫描角度的变化
%
% 函数会自动读取 PaperA_exp7_0.mat
% 数据格式：
%   - PSLR: 6×L 矩阵，6种方法的PSLR值（单位：dB）
%   - HPBW: 6×L 矩阵，6种方法的HPBW值（单位：度）
%   - scan_angle: 1×L 向量，扫描角度（单位：度）
%
% 可选参数 (通过 varargin 或结构体):
%   'figure_handle' - 图形句柄 (可选，用于在指定图形中绘制)
%   'legend_labels' - 图例标签，默认 {'Baseline', 'CDM', 'FDM', 'Random', 'Nested', 'Coprime'}
%   'line_width'   - 线宽, 默认 2.0 (科研绘图推荐)
%   'font_size'     - 字体大小, 默认 14
%   'save_fig'      - 是否保存图片, 默认 false
%   'save_path'     - 保存路径, 默认 'PaperA_exp7_PSLR_HPBW.png'
%   'fig_size'      - 图形尺寸 [width, height] (inches), 默认 [12, 5]
%   'mat_file'      - MAT文件路径, 默认 'PaperA_exp7_0.mat'
%
% 示例:
%   plot_paperA_exp7();
%   plot_paperA_exp7('save_fig', true, 'fig_size', [14, 6]);

    % ========== 参数解析 ==========
    p = inputParser;
    addParameter(p, 'figure_handle', [], @(x) isempty(x) || isgraphics(x));
    addParameter(p, 'legend_labels', {'Baseline', 'CDM', 'FDM', 'Random', 'Nested', 'Coprime'}, ...
        @(x) iscell(x) && length(x) == 6);
    addParameter(p, 'line_width', 2.0, @isnumeric);
    addParameter(p, 'font_size', 14, @isnumeric);
    addParameter(p, 'save_fig', false, @islogical);
    addParameter(p, 'save_path', 'PaperA_exp7_PSLR_HPBW.png', @(x) ischar(x) || isstring(x));
    addParameter(p, 'fig_size', [12, 5], @(x) isnumeric(x) && length(x) == 2);
    addParameter(p, 'mat_file', legacy_figure_asset('PaperA_exp7_0.mat'), @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    
    fig_handle = p.Results.figure_handle;
    legend_labels = p.Results.legend_labels;
    line_width = p.Results.line_width;
    font_size = p.Results.font_size;
    save_fig = p.Results.save_fig;
    save_path = p.Results.save_path;
    fig_size = p.Results.fig_size;
    mat_file = p.Results.mat_file;
    
    % ========== 读取MAT文件 ==========
    if ~exist(mat_file, 'file')
        error('文件不存在: %s', mat_file);
    end
    
    data = load(mat_file);
    
    % 检查必要的变量是否存在
    required_vars = {'PSLR', 'HPBW', 'scan_angle'};
    missing_vars = {};
    for i = 1:length(required_vars)
        if ~isfield(data, required_vars{i})
            missing_vars{end+1} = required_vars{i};
        end
    end
    
    if ~isempty(missing_vars)
        error('MAT文件 %s 中缺少以下变量: %s', mat_file, strjoin(missing_vars, ', '));
    end
    
    % 提取数据
    PSLR = data.PSLR;
    HPBW = data.HPBW;
    scan_angle = data.scan_angle;
    
    % 确保数据维度正确
    if size(PSLR, 1) ~= 6
        error('PSLR矩阵的行数应为6，当前为%d', size(PSLR, 1));
    end
    if size(HPBW, 1) ~= 6
        error('HPBW矩阵的行数应为6，当前为%d', size(HPBW, 1));
    end
    if size(scan_angle, 1) ~= 1 && size(scan_angle, 2) ~= 1
        error('scan_angle应为向量，当前维度为[%d, %d]', size(scan_angle, 1), size(scan_angle, 2));
    end
    
    % 确保scan_angle是行向量
    if size(scan_angle, 1) > size(scan_angle, 2)
        scan_angle = scan_angle';
    end
    
    % 检查维度一致性
    L = length(scan_angle);
    if size(PSLR, 2) ~= L
        error('PSLR矩阵的列数(%d)与scan_angle长度(%d)不匹配', size(PSLR, 2), L);
    end
    if size(HPBW, 2) ~= L
        error('HPBW矩阵的列数(%d)与scan_angle长度(%d)不匹配', size(HPBW, 2), L);
    end
    
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
    
    % ========== 创建图形 ==========
    if isempty(fig_handle)
        fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, fig_size(1), fig_size(2)]);
    else
        fig = figure(fig_handle);
    end
    
    % ========== 绘制PSLR子图 ==========
    ax1 = subplot(1, 2, 1);
    hold(ax1, 'on');
    
    h_pslr = zeros(6, 1);
    for i = 1:6
        h_pslr(i) = plot(ax1, scan_angle, PSLR(i, :), 'LineWidth', line_width, ...
            'Color', colors(i, :), 'LineStyle', line_styles{i});
    end
    
    xlabel(ax1, 'Scan Angle (degrees)', 'FontSize', font_size, 'FontWeight', 'bold');
    ylabel(ax1, 'PSLR (dB)', 'FontSize', font_size, 'FontWeight', 'bold');
    title(ax1, 'Peak Sidelobe Level Ratio', 'FontSize', font_size+1, 'FontWeight', 'bold');
    
    % 图例设置
    legend(ax1, h_pslr, legend_labels, 'Location', 'best', ...
        'FontSize', font_size-1, 'Box', 'on', 'EdgeColor', [0.3, 0.3, 0.3]);
    
    % 网格设置
    grid(ax1, 'on');
    grid(ax1, 'minor');
    ax1.GridAlpha = 0.3;
    ax1.MinorGridAlpha = 0.15;
    ax1.GridLineStyle = '-';
    ax1.MinorGridLineStyle = ':';
    
    % 坐标轴设置
    ax1.FontSize = font_size-1;
    ax1.TickLength = [0.02, 0.02];
    box(ax1, 'on');
    xlim(ax1, [min(scan_angle), max(scan_angle)]);
    
    hold(ax1, 'off');
    
    % ========== 绘制HPBW子图 ==========
    ax2 = subplot(1, 2, 2);
    hold(ax2, 'on');
    
    h_hpbw = zeros(6, 1);
    for i = 1:6
        h_hpbw(i) = plot(ax2, scan_angle, HPBW(i, :), 'LineWidth', line_width, ...
            'Color', colors(i, :), 'LineStyle', line_styles{i});
    end
    
    xlabel(ax2, 'Scan Angle (degrees)', 'FontSize', font_size, 'FontWeight', 'bold');
    ylabel(ax2, 'HPBW (degrees)', 'FontSize', font_size, 'FontWeight', 'bold');
    title(ax2, 'Half Power Beam Width', 'FontSize', font_size+1, 'FontWeight', 'bold');
    
    % 图例设置
    legend(ax2, h_hpbw, legend_labels, 'Location', 'best', ...
        'FontSize', font_size-1, 'Box', 'on', 'EdgeColor', [0.3, 0.3, 0.3]);
    
    % 网格设置
    grid(ax2, 'on');
    grid(ax2, 'minor');
    ax2.GridAlpha = 0.3;
    ax2.MinorGridAlpha = 0.15;
    ax2.GridLineStyle = '-';
    ax2.MinorGridLineStyle = ':';
    
    % 坐标轴设置
    ax2.FontSize = font_size-1;
    ax2.TickLength = [0.02, 0.02];
    box(ax2, 'on');
    xlim(ax2, [min(scan_angle), max(scan_angle)]);
    
    hold(ax2, 'off');
    
    % 调整子图间距
    set(fig, 'Units', 'inches');
    pos = get(fig, 'Position');
    set(fig, 'Position', [pos(1), pos(2), fig_size(1), fig_size(2)]);
    
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

