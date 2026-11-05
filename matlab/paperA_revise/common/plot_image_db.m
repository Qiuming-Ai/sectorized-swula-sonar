function h = plot_image_db(cfg, image_data, dynamic_range_db, plot_title)
%PLOT_IMAGE_DB Consistent angle-range dB display.
if nargin < 3 || isempty(dynamic_range_db), dynamic_range_db = 50; end
if nargin < 4, plot_title = ""; end
z = abs(image_data);
z_db = 20*log10(z/(max(z(:))+eps) + eps);
h = imagesc(cfg.angle_grid_deg, cfg.range_axis_m, z_db);
axis xy;
xlabel('Angle (deg)');
ylabel('Relative range (m)');
title(plot_title);
clim([-dynamic_range_db,0]);
colorbar;
colormap(turbo);
end
