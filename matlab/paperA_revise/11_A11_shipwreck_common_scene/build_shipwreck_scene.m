function [scene,visual_reference,ship_catalog,seabed_catalog,arrival_catalog,geometry] = build_shipwreck_scene(cfg)
%BUILD_SHIPWRECK_SCENE Build one legacy-geometry scene for all three methods.
%   The scene is derived from the one-bounce CUDA arrival cloud of
%   PaperA_exp8_baseline_real_scene.m.  The imaging coordinates also follow
%   its 1-D receive-array steering model: sin(theta)=x/r, where r is the
%   full 3-D slant range.  Using atan2(x,z) would incorrectly ignore the
%   unresolved elevation component and exaggerate the lateral angle.
%   The 45-deg direction in the old preview is not used as the imaging axis;
%   tx_cos_power=rx_cos_power=-1 disables directional gain in the ray tracer.

if ~isfile(cfg.ship_stl)
    error('A11:MissingSTL','Extract the geometry overlay or pass its geometry directory: missing %s.',cfg.ship_stl);
end

% Recreate the transformed STL for the geometry catalog and audit figure.
tr = stlread(cfg.ship_stl);
vertices_stl = cfg.ship_scale*tr.Points;
rpy = cfg.ship_rotation_rpy_deg;
rx = [1,0,0;0,cosd(rpy(1)),-sind(rpy(1));0,sind(rpy(1)),cosd(rpy(1))];
ry = [cosd(rpy(2)),0,sind(rpy(2));0,1,0;-sind(rpy(2)),0,cosd(rpy(2))];
rz = [cosd(rpy(3)),-sind(rpy(3)),0;sind(rpy(3)),cosd(rpy(3)),0;0,0,1];
rotation = rz*ry*rx;
vertices_pre_swap = (rotation*vertices_stl.').' + ...
    cfg.ship_position_before_axis_swap_m;
vertices = vertices_pre_swap(:,cfg.ship_axis_permutation);
faces = tr.ConnectivityList(1:cfg.mesh_stride:end,:);

p1 = vertices(faces(:,1),:);
p2 = vertices(faces(:,2),:);
p3 = vertices(faces(:,3),:);
cross_vec = cross(p2-p1,p3-p1,2);
double_area = sqrt(sum(cross_vec.^2,2));
centroid = (p1+p2+p3)/3;
[ship_range,ship_angle,ship_elevation,ship_forward,ship_lateral] = ...
    project_legacy_imaging_coordinates(centroid-cfg.sonar_origin_xyz_m);
ship_valid = ship_forward>0 & ship_angle>=cfg.fov_deg(1) & ...
    ship_angle<=cfg.fov_deg(2) & ship_range>=cfg.range_origin_m & ...
    ship_range<=cfg.range_origin_m+cfg.range_span_m;
centroid = centroid(ship_valid,:);
ship_range = ship_range(ship_valid);
ship_angle = ship_angle(ship_valid);
ship_elevation = ship_elevation(ship_valid);
ship_forward = ship_forward(ship_valid);
ship_lateral = ship_lateral(ship_valid);
double_area = double_area(ship_valid);
[ship_ri,ship_ai] = grid_indices(cfg,ship_range,ship_angle);
ship_catalog = table(ship_ri,ship_ai,ship_range,ship_angle,ship_elevation, ...
    centroid(:,1),centroid(:,2),centroid(:,3),0.5*double_area, ...
    ship_forward,ship_lateral, ...
    'VariableNames',{'RangeIndex','AngleIndex','RangeM','AngleDeg', ...
    'ElevationDeg','X','Y','Z','FacetArea','ForwardZM','LateralXM'});

% Recreate the exact legacy seabed point grid for the audit catalog.
rng(cfg.seabed_seed,'twister');
nx = round(sqrt(cfg.n_seabed_scatter*cfg.seabed_width_m/cfg.seabed_length_m));
nz = round(cfg.n_seabed_scatter/nx);
xlin = linspace(-cfg.seabed_width_m/2,cfg.seabed_width_m/2,nx);
zlin = linspace(0,cfg.seabed_length_m,nz);
[xg,zg] = meshgrid(xlin,zlin);
yg = cfg.seabed_y_m + cfg.seabed_roughness_sigma_m*randn(size(xg));
seabed_position = [xg(:),yg(:),zg(:)];
[seabed_range,seabed_angle,seabed_elevation,seabed_forward,seabed_lateral] = ...
    project_legacy_imaging_coordinates(seabed_position-cfg.sonar_origin_xyz_m);
seabed_valid = seabed_forward>0 & seabed_angle>=cfg.fov_deg(1) & ...
    seabed_angle<=cfg.fov_deg(2) & seabed_range>=cfg.range_origin_m & ...
    seabed_range<=cfg.range_origin_m+cfg.range_span_m;
seabed_position = seabed_position(seabed_valid,:);
seabed_range = seabed_range(seabed_valid);
seabed_angle = seabed_angle(seabed_valid);
seabed_elevation = seabed_elevation(seabed_valid);
seabed_forward = seabed_forward(seabed_valid);
seabed_lateral = seabed_lateral(seabed_valid);
[seabed_ri,seabed_ai] = grid_indices(cfg,seabed_range,seabed_angle);
seabed_catalog = table(seabed_ri,seabed_ai,seabed_range,seabed_angle, ...
    seabed_elevation,seabed_position(:,1),seabed_position(:,2), ...
    seabed_position(:,3),seabed_forward,seabed_lateral, ...
    'VariableNames',{'RangeIndex','AngleIndex','RangeM','AngleDeg', ...
    'ElevationDeg','X','Y','Z','ForwardZM','LateralXM'});

% Convert the exact one-bounce visible-arrival cloud to the common 2-D
% complex scene. Quantitative reconstruction uses the ray-tracer amplitude
% directly; display enhancement is applied only to visual_reference below.
arrival_stride = cfg.raytrace_arrival_stride;
if strcmp(cfg.arrival_source,'archive')
    if ~isfile(cfg.arrival_csv)
        error('A11:MissingArchive','Extract the geometry overlay or pass its geometry directory: missing %s.',cfg.arrival_csv);
    end
    archived = readtable(cfg.arrival_csv);
    required = {'X','Y','Z','RawAmplitude','MetricAmplitude','PathLengthM', ...
        'RangeM','AngleDeg','ElevationDeg','ForwardZM','LateralXM'};
    if ~all(ismember(required,archived.Properties.VariableNames)) || ...
            ~isequal(archived.RawAmplitude,archived.MetricAmplitude)
        error('A11:InvalidArchive','The arrival CSV does not match the common unweighted A11 arrival schema.');
    end
    % The archived CSV is already FOV/range filtered. Quick mode thins this
    % exported catalogue, not the older unfiltered CUDA arrival cache.
    archived = archived(1:arrival_stride:end,:);
    ray_position = archived{:,{'X','Y','Z'}};
    raw_amplitude = archived.RawAmplitude;
    path_length = archived.PathLengthM;
    ray_range = archived.RangeM; ray_angle = archived.AngleDeg;
    ray_elevation = archived.ElevationDeg;
    ray_forward = archived.ForwardZM; ray_lateral = archived.LateralXM;
else
    ray = audit_A11_legacy_raytrace(cfg.geometry_dir,cfg.output_dir,cfg.cuda_mex_dir);
    ray_position = double(ray.pos(1:arrival_stride:end,:));
    raw_amplitude = double(ray.amp(1:arrival_stride:end));
    path_length = double(ray.path_len(1:arrival_stride:end));
    [ray_range,ray_angle,ray_elevation,ray_forward,ray_lateral] = ...
        project_legacy_imaging_coordinates(ray_position-cfg.sonar_origin_xyz_m);
end
ray_valid = ray_forward>0 & ray_angle>=cfg.fov_deg(1) & ...
    ray_angle<=cfg.fov_deg(2) & ray_range>=cfg.range_origin_m & ...
    ray_range<=cfg.range_origin_m+cfg.range_span_m & raw_amplitude>0;
ray_position = ray_position(ray_valid,:);
raw_amplitude = raw_amplitude(ray_valid);
path_length = path_length(ray_valid);
ray_range = ray_range(ray_valid);
ray_angle = ray_angle(ray_valid);
ray_elevation = ray_elevation(ray_valid);
ray_forward = ray_forward(ray_valid);
ray_lateral = ray_lateral(ray_valid);
if isempty(ray_range)
    error('A11:EmptyROI','No positive-amplitude arrivals remain in the configured region.');
end

metric_amplitude = raw_amplitude;

% A deterministic unresolved-surface phase produces reproducible speckle.
% The separately exported visual_reference is the positive magnitude
% envelope before phase summation and is used only as a geometry reference.
rng(cfg.ship_seed,'twister');
phase = 2*pi*rand(size(metric_amplitude));
coefficient = metric_amplitude.*exp(1j*phase);
[arrival_ri,arrival_ai] = grid_indices(cfg,ray_range,ray_angle);
scene = accumarray([arrival_ri,arrival_ai],coefficient, ...
    [cfg.n_range_samples,numel(cfg.angle_grid_deg)],@sum,complex(0));
reference_envelope = accumarray([arrival_ri,arrival_ai],metric_amplitude, ...
    [cfg.n_range_samples,numel(cfg.angle_grid_deg)],@sum,0);
normalizer = max(abs(scene(:)))+eps;
scene = scene/normalizer;
sigma = cfg.reference_smoothing_sigma_bins;
range_lag = -ceil(3*sigma(1)):ceil(3*sigma(1));
angle_lag = -ceil(3*sigma(2)):ceil(3*sigma(2));
range_kernel = exp(-0.5*(range_lag/sigma(1)).^2);
angle_kernel = exp(-0.5*(angle_lag/sigma(2)).^2);
display_kernel = range_kernel(:)*angle_kernel(:).';
display_kernel = display_kernel/sum(display_kernel(:));
reference_envelope = conv2(reference_envelope,display_kernel,'same');
reference_envelope = reference_envelope/(max(reference_envelope(:))+eps);
visual_reference = reference_envelope.^cfg.reference_display_gamma;

arrival_catalog = table(arrival_ri,arrival_ai,ray_range,ray_angle, ...
    ray_elevation,raw_amplitude,metric_amplitude,phase,path_length, ...
    ray_position(:,1),ray_position(:,2),ray_position(:,3), ...
    ray_forward,ray_lateral, ...
    'VariableNames',{'RangeIndex','AngleIndex','RangeM','AngleDeg', ...
    'ElevationDeg','RawAmplitude','MetricAmplitude','PhaseRad','PathLengthM', ...
    'X','Y','Z','ForwardZM','LateralXM'});

preview_boresight = [0,sind(cfg.sonar_boresight_el_deg), ...
    cosd(cfg.sonar_boresight_el_deg)];
geometry = struct();
geometry.arrival_source = cfg.arrival_source;
geometry.preview_boresight = preview_boresight;
geometry.directionality_disabled = true;
geometry.imaging_projection = ...
    "sin(theta)=x/slant_range; slant range=sqrt(x^2+y^2+z^2)";
geometry.ship_bbox_world = [min(vertices,[],1);max(vertices,[],1)];
geometry.ship_range_limits = [min(ship_range),max(ship_range)];
geometry.ship_angle_limits = [min(ship_angle),max(ship_angle)];
geometry.ship_elevation_limits = [min(ship_elevation),max(ship_elevation)];
geometry.arrival_range_limits = [min(ray_range),max(ray_range)];
geometry.arrival_angle_limits = [min(ray_angle),max(ray_angle)];
geometry.n_ship_facets = size(ship_catalog,1);
geometry.n_seabed_points = size(seabed_catalog,1);
geometry.n_visible_ray_arrivals = size(arrival_catalog,1);
geometry.reference_display_gamma = cfg.reference_display_gamma;
geometry.reference_smoothing_sigma_bins = cfg.reference_smoothing_sigma_bins;
geometry.transform_order = ...
    "Rz(yaw)*Ry(pitch)*Rx(roll), translate, axis permutation";
end

function [range,bearing,elevation,forward,lateral] = ...
    project_legacy_imaging_coordinates(position)
range = sqrt(sum(position.^2,2));
lateral = position(:,1);
forward = position(:,3);
bearing = asind(min(max(lateral./(range+eps),-1),1));
horizontal = sqrt(lateral.^2+forward.^2);
elevation = atan2d(position(:,2),horizontal);
end

function [range_index,angle_index] = grid_indices(cfg,range,bearing)
range_index = round((range-cfg.range_origin_m)/cfg.range_sample_spacing_m)+1;
angle_step = cfg.angle_grid_deg(2)-cfg.angle_grid_deg(1);
angle_index = round((bearing-cfg.angle_grid_deg(1))/angle_step)+1;
range_index = min(max(range_index,1),cfg.n_range_samples);
angle_index = min(max(angle_index,1),numel(cfg.angle_grid_deg));
end
