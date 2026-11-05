% TEST_PLOT_EXP7 测试 plot_paperA_exp7 函数
% 快速测试脚本，用于验证可视化功能

clear; close all; clc;

% 检查文件是否存在
mat_file = legacy_figure_asset('PaperA_exp7_0.mat');
if ~exist(mat_file, 'file')
    error('文件不存在: %s', mat_file);
end

% 读取数据并显示基本信息
data = load(mat_file);
fprintf('文件 %s 中的数据变量:\n', mat_file);
fprintf('  PSLR size: [%d, %d]\n', size(data.PSLR));
fprintf('  HPBW size: [%d, %d]\n', size(data.HPBW));
fprintf('  scan_angle size: [%d, %d]\n', size(data.scan_angle));
fprintf('\n');

% 调用可视化函数
fprintf('正在绘制图形...\n');
plot_paperA_exp7('save_fig', false);

fprintf('可视化完成！\n');

