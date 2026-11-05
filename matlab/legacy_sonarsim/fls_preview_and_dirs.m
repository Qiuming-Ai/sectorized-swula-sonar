function [tx_dir, rx_dir] = fls_preview_and_dirs(models, S, R, opts)
% FLS 场景预览 + 方向生成 + 可视化
% 新增：
%   opts.model_color_mode : 'byModel' | 'byReflect'  （默认 'byModel'：不同模型不同色）
%   opts.model_colors     : [Kx3] 自定义RGB调色板（0~1），少于模型数则循环使用
%   opts.face_alpha       : 面透明度（默认 0.35）
%   opts.edge_alpha       : 边线透明度（默认 0.35）
%   opts.legend           : 是否显示图例（默认 true）
%
% 使用示例见文件末尾 demo()

% ------------------ 参数 ------------------
if nargin<4 || isempty(opts) || ~isstruct(opts), opts = struct; end
opts = setdef(opts, 'tx_az_deg', 0);           % 发射方位 (deg)
opts = setdef(opts, 'tx_el_deg', -5);          % 发射俯仰 (deg), 例如 -5° 俯视
opts = setdef(opts, 'rx_mode', 'common');      % 'common' | 'toFocus' | 'sector'
opts = setdef(opts, 'rx_el_deg', 0);           % rx 公共俯仰（common/sector 用）
opts = setdef(opts, 'rx_common_az_deg', 0);    % rx 公共方位（common）
opts = setdef(opts, 'rx_focus', []);           % [1x3] 焦点（toFocus 用）
opts = setdef(opts, 'rx_sector_span_deg', 80); % 扇扫总宽（sector 用）
opts = setdef(opts, 'rx_sector_center_deg', 0);% 扇扫中心方位（sector 用）
opts = setdef(opts, 'draw_rx_arrows', true);
opts = setdef(opts, 'draw_hemisphere', true);  % 画半球门控示意
opts = setdef(opts, 'hemi_radius', []);
opts = setdef(opts, 'hemi_density', 18);       % 半球网格密度
opts = setdef(opts, 'fig', []);
opts = setdef(opts, 'axis_equal', true);
opts = setdef(opts, 'bg_color', [1 1 1]);

% —— 新增可选项 —— 
opts = setdef(opts, 'model_color_mode', 'byModel'); % 'byModel' | 'byReflect'
opts = setdef(opts, 'model_colors', []);            % [Kx3], 0~1
opts = setdef(opts, 'face_alpha', 0.35);
opts = setdef(opts, 'edge_alpha', 0.35);
opts = setdef(opts, 'legend', true);

Nr = size(R,1);

% ------------------ 方向生成 ------------------
tx_dir = dir_from_az_el(opts.tx_az_deg, opts.tx_el_deg);

switch lower(opts.rx_mode)
    case 'common'
        dirc = dir_from_az_el(opts.rx_common_az_deg, opts.rx_el_deg);
        rx_dir = repmat(dirc, Nr, 1);
    case 'tofocus'
        if isempty(opts.rx_focus), error('rx_mode=toFocus 需要 opts.rx_focus=[1x3]'); end
        F = opts.rx_focus(:).';
        v  = F - R;  rx_dir = normalize_rows(v);
    case 'sector'
        span = opts.rx_sector_span_deg;
        cen  = opts.rx_sector_center_deg;
        azs  = linspace(cen - span/2, cen + span/2, Nr);
        rx_dir = zeros(Nr,3);
        for i=1:Nr
            rx_dir(i,:) = dir_from_az_el(azs(i), opts.rx_el_deg);
        end
    otherwise
        error('未知 rx_mode: %s', opts.rx_mode);
end

% 半球半径默认：取场景尺度的 15%
if isempty(opts.hemi_radius)
    box = scene_bbox(models, S, R);
    diagL = norm(box(2,:) - box(1,:));
    opts.hemi_radius = 0.15*diagL;
end

% ------------------ 可视化 ------------------
if isempty(opts.fig) || ~isvalid(opts.fig)
    f = figure('Color', opts.bg_color, 'Name','FLS Scene Preview');
else
    f = opts.fig; figure(f); clf(f);
end
ax = axes('Parent', f); hold(ax, 'on');

% —— 生成调色板（按模型）——
Nmod = numel(models);
palette = [];
if strcmpi(opts.model_color_mode,'byModel')
    if isempty(opts.model_colors)
        palette = lines(max(Nmod,7)); % 稳健、对比友好的默认调色板
    else
        palette = opts.model_colors;
        if size(palette,2)~=3, error('opts.model_colors 需为 [Kx3] RGB（0~1）'); end
    end
end

% 画 STL 模型（不同模型不同颜色 / 或按反射率着色）
ph = gobjects(Nmod,1);
legtxt = strings(Nmod,1);
for m = 1:Nmod
    M = models{m};
    if ~isfield(M,'Points') || ~isfield(M,'ConnectivityList') || isempty(M.Points)
        warning('models{%d} 缺少 Points/ConnectivityList 或为空，已跳过。', m);
        continue;
    end
    V = M.Points; 
    Fcs = M.ConnectivityList;

    switch lower(opts.model_color_mode)
        case 'bymodel'
            c = palette( mod(m-1, size(palette,1)) + 1, : );
            fcol = c;
            ecol = 0.55*c + 0.25; % 略亮的边色
        case 'byreflect'
            rf = 0.5;
            if isfield(M,'reflectivity') && ~isempty(M.reflectivity)
                rf = min(max(M.reflectivity,0),1);
            end
            % 用 colormap (parula) 取色：低反射→蓝，高反射→黄
            cmap = parula(256);
            idx = 1 + round(rf * 255);
            fcol = cmap(idx,:);
            ecol = 0.6*fcol + 0.25;
        otherwise
            error('未知 opts.model_color_mode: %s', opts.model_color_mode);
    end

    ph(m) = patch(ax, 'Faces', Fcs, 'Vertices', V, ...
        'FaceColor', fcol, 'EdgeColor', ecol, ...
        'FaceAlpha', opts.face_alpha, 'EdgeAlpha', opts.edge_alpha);

    if isfield(M,'name') && ~isempty(M.name)
        legtxt(m) = string(M.name);
    else
        legtxt(m) = "Model " + m;
    end
end

% 画 S、R
plot3(ax, S(1),S(2),S(3), 'rp', 'MarkerFaceColor','r', 'MarkerSize',10);
plot3(ax, R(:,1),R(:,2),R(:,3), 'bo', 'MarkerFaceColor',[0.2 0.5 1], 'MarkerSize',5);

% 发射箭头
qscale = opts.hemi_radius*0.9;
quiver3(ax, S(1),S(2),S(3), tx_dir(1)*qscale, tx_dir(2)*qscale, tx_dir(3)*qscale, ...
    0, 'LineWidth',2.0, 'Color',[0.85 0.2 0.2]);

% 接收箭头
if opts.draw_rx_arrows
    U = rx_dir * qscale;
    quiver3(ax, R(:,1),R(:,2),R(:,3), U(:,1),U(:,2),U(:,3), 0, ...
        'LineWidth', 1.2, 'Color',[0.2 0.4 0.9]);
end

% 半球门控示意（发射端 + 可选一个接收端）
if opts.draw_hemisphere
    draw_hemi(ax, S, tx_dir, opts.hemi_radius, opts.hemi_density, [1 0.6 0.6], 0.15);
    if size(R,1)>0
        draw_hemi(ax, R(1,:), rx_dir(1,:), opts.hemi_radius*0.7, opts.hemi_density, [0.6 0.7 1], 0.15);
    end
end

% 轴、视图 & 图例
axis(ax, 'tight'); grid(ax, 'on'); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
if opts.axis_equal, axis(ax,'equal'); end
view(ax, 45, 25);
title(ax, sprintf('tx az=%.1f°, el=%.1f°  |  rx: %s', ...
    opts.tx_az_deg, opts.tx_el_deg, upper(opts.rx_mode)));

if opts.legend
    valid = isgraphics(ph);
    if any(valid)
        legend(ax, ph(valid), legtxt(valid), 'Location','northeastoutside', ...
            'Interpreter','none', 'AutoUpdate','off');
    end
end

end % main

% ===================== 辅助函数 =====================
function s = setdef(s, k, defaultVal)
    if ~isstruct(s), s = struct; end
    if ~isfield(s, k) || isempty(s.(k))
        s.(k) = defaultVal;
    end
end

function d = dir_from_az_el(az_deg, el_deg)
% 坐标系：+z 前方，+x 右，+y 上
    az = deg2rad(az_deg);
    el = deg2rad(el_deg);
    d = [sin(az).*cos(el),  sin(el),  cos(az).*cos(el)];
    d = d ./ vecnorm(d,2,2);
end

function Vn = normalize_rows(V)
    n = sqrt(sum(V.^2,2))+eps;
    Vn = V ./ n;
end

function bbox = scene_bbox(models, S, R)
    mins = [inf inf inf];  maxs = -mins;
    for m=1:numel(models)
        if isfield(models{m},'Points') && ~isempty(models{m}.Points)
            V = models{m}.Points;
            mins = min(mins, min(V,[],1));
            maxs = max(maxs, max(V,[],1));
        end
    end
    if ~isempty(S), mins = min(mins, S); maxs = max(maxs, S); end
    if ~isempty(R), mins = min(mins, min(R,[],1)); maxs = max(maxs, max(R,[],1)); end
    bbox = [mins; maxs];
end

function draw_hemi(ax, C, dir, r, dens, color, alphaV)
% 画以 dir 为轴向、半径 r 的半球网格（稳健版）
    d = dir(:).';
    nd = norm(d);
    if nd < 1e-12, d = [0 0 1]; else, d = d / nd; end

    UV = null(d).';
    if size(UV,1) >= 2
        u = UV(1,:); v = UV(2,:);
    else
        if abs(d(1)) < 0.9, tmp = [1 0 0]; else, tmp = [0 1 0]; end
        u = tmp - dot(tmp,d)*d;     u = u / norm(u);
        v = cross(d,u);             v = v / norm(v);
    end
    u = u / norm(u);
    v = v - dot(v,u)*u; v = v / norm(v);

    th = linspace(0, pi/2, dens);
    ph = linspace(0, 2*pi, max(2*dens,12));
    [PH, TH] = meshgrid(ph, th);
    Xl = sin(TH).*cos(PH);  Yl = sin(TH).*sin(PH);  Zl = cos(TH);

    X = C(1) + r*( Xl*u(1) + Yl*v(1) + Zl*d(1) );
    Y = C(2) + r*( Xl*u(2) + Yl*v(2) + Zl*d(2) );
    Z = C(3) + r*( Xl*u(3) + Yl*v(3) + Zl*d(3) );

    surf(ax, X,Y,Z, 'FaceColor', color, 'FaceAlpha', alphaV, ...
         'EdgeAlpha', 0.2, 'EdgeColor',[0.5 0.5 0.5]);
end

% ===================== 演示 =====================
function demo()
% 假设已有 models, S, R 变量：
% (1) 不同模型不同颜色（默认）：
% [tx_dir, rx_dir] = fls_preview_and_dirs(models, S, R, struct( ...
%   'tx_az_deg', 0, 'tx_el_deg', -6, ...
%   'rx_mode','sector', 'rx_sector_span_deg', 90, 'rx_sector_center_deg', 0, ...
%   'rx_el_deg', 0, 'legend', true));
%
% (2) 自定义调色板（RGB）：
% mycolors = [0.90 0.33 0.32; 0.49 0.78 0.28; 0.30 0.55 0.90; 0.93 0.69 0.13];
% [tx_dir, rx_dir] = fls_preview_and_dirs(models, S, R, struct('model_colors', mycolors));
%
% (3) 按反射率着色：
% [~,~] = fls_preview_and_dirs(models, S, R, struct('model_color_mode','byReflect'));
end
