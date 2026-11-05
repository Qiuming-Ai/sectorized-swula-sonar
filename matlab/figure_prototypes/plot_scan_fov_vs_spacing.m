function plot_scan_fov_vs_spacing()
% Bands only: FOV (±theta_max), mainlobe width (±0.5 HPBW at theta0),
% and grating-lobe widths (m=1,2) computed at their true centers
% sin(theta_g) = sin(theta0) ± m/x, visible when |sin(theta_g)| ≤ 1.
%
% x-axis: d/lambda in [0, 2.5]
% Text is parser-safe (Interpreter='none'). No lines are drawn.

% ---------------- User params ----------------
theta0_deg = 0;     % steering angle theta0 in degrees (change if needed)
N_elt      = 16;    % number of elements (affects HPBW)
x_max      = 2.5;   % d/lambda upper bound
orders     = [1 2]; % grating orders to show

% ---------------- Internal setup ----------------
x_min_calc = 0.02;    % start >0 to avoid division by zero
Nx          = 2200;   % sampling
x = linspace(x_min_calc, x_max, Nx);
s0 = sind(theta0_deg);

% Colors
col_fov  = [0.75 0.85 1.00]; % light blue
col_main = [0.20 0.60 0.20]; % green
col_m1   = [0.90 0.55 0.15]; % amber
col_m2   = [0.55 0.35 0.95]; % violet

% Constant for HPBW (deg): 0.886 * (lambda / (L * cos)) * 180/pi, L ~ (N-1)*d
K_deg = 0.886 * (180/pi);

% ------------- FOV: theta_max vs x -------------
theta_max = zeros(size(x));
mask1 = (x <= 0.5);
mask2 = (x > 0.5 & x < 1.0);
mask3 = (x >= 1.0);
theta_max(mask1) = 90;
theta_max(mask2) = asind( max(0, 1./x(mask2) - 1) );
theta_max(mask3) = 0;

% ------------- Mainlobe width at theta0 -------------
HPBW_main_deg = K_deg ./ ((N_elt-1) .* x .* max(cosd(theta0_deg), 0.02));
half_main = min(0.5 * HPBW_main_deg, 90);

% ------------- Figure (bands only) -------------
figure('Color','w'); hold on;

% FOV band (±theta_max)
Xfill = [x, fliplr(x)];
Yfov  = [theta_max, -fliplr(theta_max)];
pFOV = fill(Xfill, Yfov, col_fov, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
    'DisplayName', 'FOV (grating-free, +/-theta_max)');

% Mainlobe width band centered at theta0
Ymain = [theta0_deg + half_main, fliplr(theta0_deg - half_main)];
pMain = fill(Xfill, Ymain, col_main, 'FaceAlpha', 0.28, 'EdgeColor', 'none', ...
    'DisplayName', sprintf('Mainlobe width at %.1f deg (+/-0.5 HPBW)', theta0_deg));

% ------------- Grating-lobe width bands (m=1,2) -------------
pGL = gobjects(0);
for m = orders
    for sgn = [+1, -1]
        % Center sine: sin(theta_g) = sin(theta0) + sgn * m / x
        u = s0 + sgn * (m ./ x);         % may be outside [-1,1]
        vis = (abs(u) <= 1);             % visible region for this sign/order
        if ~any(vis), continue; end

        x_vis = x(vis);
        theta_g = asind(u(vis));         % center angle (deg)
        cos_g   = max(cosd(theta_g), 0.02);

        % HPBW at the grating-lobe center
        HPBW_gl_deg = K_deg ./ ((N_elt-1) .* x_vis .* cos_g);
        half_gl = min(0.5 * HPBW_gl_deg, 90);

        Ytop = min( 90, theta_g + half_gl);
        Ybot = max(-90, theta_g - half_gl);

        % color per order
        col = col_m1; if m==2, col = col_m2; end

        % fill band for this sign/order
        h = fill([x_vis, fliplr(x_vis)], [Ytop, fliplr(Ybot)], col, ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
            'DisplayName', sprintf('Grating-lobe width m=%d (+/-0.5 HPBW)', m));

        % only keep one legend entry per order
        if sgn == -1
            set(h, 'HandleVisibility', 'off');
        else
            pGL(end+1) = h; %#ok<AGROW>
        end
    end
end

% ------------- Axes & legend -------------
xlim([0, x_max]); ylim([-90, 90]);
grid on; box on;
xlabel('Element spacing / Wavelength (d/lambda)', 'Interpreter','none');
ylabel('Angle (degrees)', 'Interpreter','none');
title(sprintf('Bands: FOV, Mainlobe, and Grating-lobe Widths vs d/lambda (N=%d, theta0=%.1f deg)', ...
      N_elt, theta0_deg), 'Interpreter','none');

if ~isempty(pGL)
    legend([pFOV, pMain, pGL(:)'], 'Location','northoutside', ...
           'NumColumns', 3, 'Interpreter','none');
else
    legend([pFOV, pMain], 'Location','northoutside', ...
           'NumColumns', 2, 'Interpreter','none');
end

% ------------- Save -------------
outfile = 'bands_FOV_main_GL_corrected_0_to_2p5lambda.png';
exportgraphics(gcf, outfile, 'Resolution', 200);
fprintf('Saved: %s\n', outfile);
end
