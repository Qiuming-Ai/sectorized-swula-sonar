function cfg = config_A11(mode, geometry_dir, arrival_source)
%CONFIG_A11 Common-scene shipwreck channel experiment in the legacy geometry.
%   The ship pose, scale, position, axis permutation, seabed, CUDA-visible
%   arrivals, and x-z receive-array imaging plane follow the baseline Exp-8
%   script. Dense/CDM/FDM share one projected scene and processing chain.
if nargin<1 || isempty(mode), mode="quick"; end
mode = validatestring(mode,{'quick','paper'});
root = fileparts(fileparts(mfilename('fullpath')));
matlab_root = fileparts(root);
package_root = fileparts(matlab_root);
if nargin<2 || isempty(geometry_dir)
    geometry_dir = fullfile(package_root,'data','shipwreck','geometry');
end
if nargin<3 || isempty(arrival_source), arrival_source="archive"; end
arrival_source = validatestring(arrival_source,{'archive','retrace'});
addpath(fullfile(root,'common'));
cfg = exp_default_config(mode);
cfg.methods = ["dense","cdm","fdm_current"];
cfg.method_labels = ["Dense ULA","CDM SW-ULA","FDM SW-ULA"];
cfg.n_rx_dense = 128;
cfg.waveform_fs = 200e3;
cfg.range_sample_spacing_m = cfg.c/(2*cfg.waveform_fs);
% The legacy ship is centred at [0,40,40] m and spans approximately
% 46.4--117.1 m in slant range after its original transform.  The 4.5--6 m
% plotting range in the legacy script is inconsistent with this geometry.
cfg.range_origin_m = 42.0;
cfg.range_span_m = 80.0;
cfg.n_range_samples = round(cfg.range_span_m/cfg.range_sample_spacing_m)+1;
cfg.range_axis_m = (0:cfg.n_range_samples-1)*cfg.range_sample_spacing_m;
cfg.physical_range_axis_m = cfg.range_origin_m+cfg.range_axis_m;
if lower(string(mode))=="quick"
    cfg.angle_grid_deg = cfg.fov_deg(1):1:cfg.fov_deg(2);
    cfg.mesh_stride = 8;
    cfg.n_seabed_scatter = 10000;
    cfg.raytrace_arrival_stride = 8;
else
    cfg.angle_grid_deg = cfg.fov_deg(1):0.5:cfg.fov_deg(2);
    cfg.mesh_stride = 1;
    cfg.n_seabed_scatter = 100000;
    cfg.raytrace_arrival_stride = 1;
end
cfg.geometry_dir = geometry_dir;
cfg.ship_stl = fullfile(geometry_dir,'ship.stl');
cfg.arrival_csv = fullfile(geometry_dir,'A11_raytrace_arrivals.csv');
cfg.arrival_source = arrival_source;
cfg.cuda_mex_dir = fullfile(package_root,'outputs','matlab','shipwreck_cuda');
cfg.output_dir = fullfile(package_root,'outputs','matlab',"A11_"+string(mode)+"_"+string(arrival_source));
cfg.ship_scale = 1.0;
cfg.ship_position_before_axis_swap_m = [0,40,40];
cfg.ship_rotation_rpy_deg = [90,0,-45];
cfg.ship_axis_permutation = [1,3,2];
cfg.sonar_origin_xyz_m = [0,0,0];
cfg.sonar_boresight_az_deg = 0;
cfg.sonar_boresight_el_deg = 45;
cfg.ship_reflectivity = 0.7;
cfg.ship_scattering_coeff = 1.0;
cfg.ship_seed = 20250821;
cfg.reference_display_gamma = 0.45;
cfg.reference_smoothing_sigma_bins = [1.0,0.35];

% Legacy seabed: x width 200 m, forward z length 100 m, depth/elevation
% coordinate y=48 m, with 0.1-m Gaussian roughness.
cfg.seabed_width_m = 200;
cfg.seabed_length_m = 100;
cfg.seabed_y_m = 48;
cfg.seabed_roughness_sigma_m = 0.1;
cfg.seabed_reflectivity = 0.1;
cfg.seabed_scattering_coeff = 0.9;
cfg.seabed_seed = 20250822;
cfg.input_snr_db = 10;
cfg.noise_seed = 2026081402;
cfg.metric_guard_radius = [16,2];
cfg.description = "Historical MATLAB shipwreck control; not the final Python processing or statistics in the revised manuscript";
end
