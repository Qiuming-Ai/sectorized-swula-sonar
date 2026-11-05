function out = audit_A11_legacy_raytrace(geometry_dir,out_dir,mex_dir)
%AUDIT_A11_LEGACY_RAYTRACE Reproduce the Exp-8 one-bounce arrival cloud.
%   This helper keeps the original ship transform, seabed triangulation,
%   material parameters, and CUDA ray tracer.  It deliberately uses one
%   co-located receiver because the output is a common scene description,
%   not method-specific channel data.

if nargin<1, geometry_dir=[]; end
cfg = config_A11("paper",geometry_dir,"retrace");
if nargin<2 || isempty(out_dir), out_dir=cfg.output_dir; end
if nargin<3 || isempty(mex_dir), mex_dir=cfg.cuda_mex_dir; end
mex_name = 'stl_raytrace_arrivals_cuda3_5_1';
if ~isfile(fullfile(mex_dir,[mex_name '.' mexext]))
    error('A11:MissingMEX','Build the local MEX with matlab/shipwreck_cuda/build_shipwreck_mex.m first; expected %s.',mex_dir);
end
if ~isfile(cfg.ship_stl)
    error('A11:MissingSTL','Extract the geometry overlay or pass its geometry directory: missing %s.',cfg.ship_stl);
end
old_path = path;
cleanup_path = onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(mex_dir,'-begin');

c0 = 1500;
f0 = 0.45e6;
lambda = c0/f0;
origin = [0,0,0];

ship_tr = stlread(cfg.ship_stl);
ship_vertices = ship_tr.Points;
r = [90,0,-45];
rx = [1,0,0;0,cosd(r(1)),-sind(r(1));0,sind(r(1)),cosd(r(1))];
ry = [cosd(r(2)),0,sind(r(2));0,1,0;-sind(r(2)),0,cosd(r(2))];
rz = [cosd(r(3)),-sind(r(3)),0;sind(r(3)),cosd(r(3)),0;0,0,1];
ship_vertices = (rz*ry*rx*ship_vertices.').' + [0,40,40];
ship_vertices = ship_vertices(:,[1,3,2]);

n_seabed = 1e5;
seabed_len = 100;
seabed_wid = 200;
nx = round(sqrt(n_seabed*seabed_wid/seabed_len));
nz = round(n_seabed/nx);
xlin = linspace(-seabed_wid/2,seabed_wid/2,nx);
zlin = linspace(0,seabed_len,nz);
[xg,zg] = meshgrid(xlin,zlin);
rng(20250822,'twister');
yg = 48 + 0.1*randn(size(xg));
seabed_vertices = [xg(:),yg(:),zg(:)];
[nz,nx] = size(xg);
index = reshape(1:nz*nx,nz,nx);
v1 = index(1:end-1,1:end-1); v2 = index(1:end-1,2:end);
v3 = index(2:end,1:end-1);   v4 = index(2:end,2:end);
seabed_faces = [[v1(:),v2(:),v4(:)];[v1(:),v4(:),v3(:)]];

ship_model = struct( ...
    'Points',ship_vertices, ...
    'ConnectivityList',ship_tr.ConnectivityList, ...
    'reflectivity',0.7, ...
    'scattering_coeff',1, ...
    'center_frequency',f0, ...
    'transmittance',0.1, ...
    'isHollow',false, ...
    'internal_scattering',0.5);
seabed_model = struct( ...
    'Points',seabed_vertices, ...
    'ConnectivityList',seabed_faces, ...
    'reflectivity',0.1, ...
    'scattering_coeff',0.9, ...
    'center_frequency',f0, ...
    'transmittance',0.8, ...
    'isHollow',true, ...
    'internal_scattering',0);

look = [0,sind(45),cosd(45)];
params = struct();
params.f0 = f0;
params.c = c0;
params.rx_dir = look;
params.tx_dir = look;
params.max_bounces = 1;
params.samples_per_tri = 2;
params.eps = 1e-5;
params.seed = uint64(20250821);
params.tx_cos_power = -1;
params.rx_cos_power = -1;
params.brdf_cos_power = 1.0;
params.rough_sigma = lambda*16;
params.amp_type = 0;
params.reflect_is_power = 1;
params.order2_extra = 0.7;
params.order3_extra = 0.4;
params.max_arrivals = 1e6;

out = stl_raytrace_arrivals_cuda3_5_1( ...
    {ship_model,seabed_model},origin,origin,params);
out.ship_vertices = ship_vertices;
out.ship_faces = ship_tr.ConnectivityList;
out.seabed_vertices = seabed_vertices;
out.seabed_faces = seabed_faces;
out.slant_range_m = sqrt(sum(out.pos.^2,2));
out.imaging_angle_deg = asind(min(max(out.pos(:,1)./(out.slant_range_m+eps),-1),1));

if ~exist(out_dir,'dir'), mkdir(out_dir); end
save(fullfile(out_dir,'A11_legacy_raytrace_cache.mat'),'out','-v7.3');
fprintf('Legacy ray trace: %d arrivals, range %.3f--%.3f m, angle %.3f--%.3f deg.\n', ...
    size(out.pos,1),min(out.slant_range_m),max(out.slant_range_m), ...
    min(out.imaging_angle_deg),max(out.imaging_angle_deg));
end
