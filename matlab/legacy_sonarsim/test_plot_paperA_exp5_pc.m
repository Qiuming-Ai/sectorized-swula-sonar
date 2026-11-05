% 测试绘制 PaperA_exp5 脉冲压缩结果（科研绘图标准）
% 函数会自动读取 PaperA_exp5_0_pc.mat
clear all; close all;
addpath(genpath("./Lib"))

% 基本调用（叠加模式：所有三条线在同一图中）
fprintf('绘制 PaperA_exp5 脉冲压缩结果（科研绘图标准）...\n');
fprintf('读取文件: PaperA_exp5_0_pc.mat\n');
plot_paperA_exp5_pulsecompression();

% 可选：使用子图模式（每个编码方式单独一个子图）
% figure;
% plot_paperA_exp5_pulsecompression('plot_type', 'subplot');

% 可选：使用自定义参数并保存图片
% figure;
% plot_paperA_exp5_pulsecompression( ...
%     'dyn_range', 40, ...
%     'line_width', 2.5, ...
%     'font_size', 16, ...
%     'save_fig', true, ...
%     'save_path', 'pulsecompression_paperA_exp5.png', ...
%     'fig_size', [12, 8], ...
%     'plot_type', 'overlay');

fprintf('完成！\n');










