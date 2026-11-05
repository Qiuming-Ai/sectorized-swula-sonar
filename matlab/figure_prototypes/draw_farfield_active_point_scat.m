function draw_farfield_active_point_scat(tx_x, rx_x, x0_tx, x0_rx, S, m_idx, n_idx, grid_center)
% 远场主动点散射示意（阵列版）
% 坐标轴：角度-距离 (θ–r)，0° 取原 y 轴正半轴方向，角度顺时针为正
% 只扫描：-75° ~ +75° 扇区；网格加密
%
% 符号标注：
%   发射/接收阵元：x_m^{tx}=(x_m^{tx},0), x_n^{rx}=(x_n^{rx},0)
%   参考中心：x_0^{tx}=(x_0^{tx},0), x_0^{rx}=(x_0^{rx},0)
%
% 输入（均可省略，使用默认）：
%   tx_x [1×M] 发射阵元 x 坐标；rx_x [1×N] 接收阵元 x 坐标
%   x0_tx, x0_rx 参考中心 x 坐标
%   S=[x_s,y_s]  散射点
%   m_idx, n_idx 高亮的发射/接收阵元索引
%   grid_center=[xc,yc] 极点（默认 (x0_tx+x0_rx)/2, 0）

%% -------- 默认参数 --------
if nargin < 1 || isempty(tx_x),  tx_x  = linspace(-0.8, 0.8, 16); end
if nargin < 2 || isempty(rx_x),  rx_x  = linspace( 1.2, 2.8, 16); end
if nargin < 3 || isempty(x0_tx), x0_tx = mean(tx_x);              end
if nargin < 4 || isempty(x0_rx), x0_rx = mean(rx_x);              end
if nargin < 5 || isempty(S),     S     = [1.0, 1.1];              end
if nargin < 6 || isempty(m_idx), m_idx = ceil(numel(tx_x)/2);     end
if nargin < 7 || isempty(n_idx), n_idx = ceil(numel(rx_x)/2);     end
if nargin < 8 || isempty(grid_center)
    grid_center = [ (x0_tx + x0_rx)/2 , 0 ];                      % 极点
end
figure
%% -------- 颜色 --------
txCol   = [0.00 0.45 0.74];   % Tx 蓝
rxCol   = [0.85 0.33 0.10];   % Rx 橙
ctrCol  = [0.25 0.25 0.25];   % 参考中心 灰
scatCol = [0.00 0.60 0.50];   % 散射点 青
rayCol  = [0.15 0.15 0.15];   % 路径 深灰

%% -------- 数据准备（极坐标转换） --------
Tx = [tx_x(:), zeros(numel(tx_x),1)];
Rx = [rx_x(:), zeros(numel(rx_x),1)];
F1 = Tx(m_idx,:);                      % 选中的 Tx 元素
F2 = Rx(n_idx,:);                      % 选中的 Rx 元素

% 以 grid_center 为极点，将 (x,y) -> (theta,r)
to_polar = @(P) cart2pol(P(:,1)-grid_center(1), P(:,2)-grid_center(2)); % theta(rad) from +x, CCW
% 我们显示时设置：ThetaZeroLocation='top'(0°=+y)，ThetaDir='clockwise'(顺时针为正)

[th_tx, r_tx] = to_polar(Tx);
[th_rx, r_rx] = to_polar(Rx);
[th_F1, r_F1] = to_polar(F1);
[th_F2, r_F2] = to_polar(F2);
[th_S,  r_S ] = to_polar(S);

% 半径范围
r_all = [r_tx; r_rx; r_S];
r_max = max(r_all)*1.15 + 1e-6;

%% -------- 极坐标轴（θ–r） --------
ax = polaraxes('Color','w');
ax.ThetaZeroLocation = 'top';        % 0° 在 +y
ax.ThetaDir          = 'clockwise';  % 顺时针为正
ax.ThetaLim          = [-75 75];     % 只扫描 -75°~75°
ax.RLim              = [0 r_max];

% 加密网格：角向 5°，半径 8-10 圈
ax.ThetaTick = -75:5:75;
nr = 9;                               % 同心圈数量-1（可调）
ax.RTick = linspace(0, r_max, nr+1);

% 如需隐藏数值标签（保持纯符号示意），取消注释下一行：
ax.ThetaTickLabel = []; ax.RTickLabel = [];

ax.GridAlpha = 0.35;      % 网格透明度
ax.GridColor = [0.8 0.8 0.8];

hold(ax,'on');

%% -------- 绘制阵列元素（θ–r） --------
% 注意：若某些阵元位于 ±90°，会超出 ±75° 的显示范围而被裁剪，这是预期的
h_tx = polarplot(ax, th_tx, r_tx, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor',txCol,'MarkerFaceColor',txCol*0.1+0.9);
h_rx = polarplot(ax, th_rx, r_rx, 's', ...
    'MarkerSize',5,'MarkerEdgeColor',rxCol,'MarkerFaceColor',rxCol*0.1+0.9);

% 高亮第 m、n 个阵元
polarplot(ax, th_F1, r_F1, 'o', 'MarkerSize',9, 'MarkerFaceColor',txCol,'MarkerEdgeColor',txCol);
polarplot(ax, th_F2, r_F2, 's', 'MarkerSize',9, 'MarkerFaceColor',rxCol,'MarkerEdgeColor',rxCol);

% 参考中心（以极点显示）
text(ax, 0, 0, '$x_0^{\mathrm{tx}},\,x_0^{\mathrm{rx}}$', 'Interpreter','latex', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Color',ctrCol);

%% -------- 散射点与路径（采样直线后映射到极坐标绘制） --------
polarplot(ax, th_S, r_S, 'p', 'MarkerSize',11, 'MarkerFaceColor',scatCol,'MarkerEdgeColor',scatCol);

% 直线 F1 -> S, S -> F2 以笛卡尔采样再映射，确保显示为直观“射线”
plot_cart_line_in_polar(ax, F1, S,  grid_center, rayCol, 180);
plot_cart_line_in_polar(ax, S,  F2, grid_center, rayCol, 180);

%% -------- 符号标注（θ–r 坐标系） --------
% 文字位置做一点径向偏移，保持清晰
rofs = 0.04*r_max;
text(ax, th_F1, r_F1+rofs, '$x_m^{\mathrm{tx}}$', 'Interpreter','latex', ...
    'HorizontalAlignment','center','Color',txCol,'FontWeight','bold');
text(ax, th_F2, r_F2+rofs, '$x_n^{\mathrm{rx}}$', 'Interpreter','latex', ...
    'HorizontalAlignment','center','Color',rxCol,'FontWeight','bold');
text(ax, th_S,  r_S +rofs, '散射点', 'HorizontalAlignment','center', 'Color',scatCol,'FontWeight','bold');

title(ax, '远场主动点散射示意（阵列版 · θ–r 扇区：-75^\circ\sim 75^\circ）', 'Interpreter','tex');

legend(ax, [h_tx, h_rx], {'Tx 阵元','Rx 阵元'}, ...
    'Location','southoutside','NumColumns',2,'Box','off');

end

%% ====== 工具：在极坐标轴上绘制笛卡尔直线（采样后映射） ======
function plot_cart_line_in_polar(ax, P1, P2, center, col, N)
if nargin<6 || isempty(N), N = 200; end
t  = linspace(0,1,N).';
XY = P1 + t.*(P2 - P1);
[th, r] = cart2pol(XY(:,1)-center(1), XY(:,2)-center(2));
polarplot(ax, th, r, '-', 'Color', col, 'LineWidth', 1.8);
end
