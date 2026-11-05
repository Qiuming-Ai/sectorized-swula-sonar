% 测试绘制 PaperA_exp5 波束图（科研绘图标准）
% 函数会自动读取 PaperA_exp5_0.mat 和 PaperA_exp5_1.mat
clear all; close all;
addpath(genpath("./Lib"))

% 基本调用（使用默认参数：Baseline, CDM, FDM, Random, Nested, Coprime）
fprintf('绘制 PaperA_exp5 波束图（科研绘图标准）...\n');
fprintf('读取文件: PaperA_exp5_0.mat 和 PaperA_exp5_1.mat\n');
plot_paperA_exp5_beampattern();

% 可选：使用自定义参数并保存图片
% figure;
% plot_paperA_exp5_beampattern( ...
%     'dyn_range', 40, ...
%     'line_width', 2.5, ...
%     'font_size', 16, ...
%     'save_fig', true, ...
%     'save_path', 'beampattern_paperA_exp5.png', ...
%     'fig_size', [12, 8]);

fprintf('完成！\n');

