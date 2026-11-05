function [scene, catalog] = generate_scatter_scene(cfg, n_scatter, n_layers, dynamic_range_db, seed)
%GENERATE_SCATTER_SCENE Reproducible sparse/dense complex reflectivity map.
if nargin < 5
    seed = cfg.default_seed;
end
rng(seed, 'twister');
nr = cfg.n_range_samples;
na = numel(cfg.angle_grid_deg);
n_layers = min(max(1, round(n_layers)), nr-8);

layer_idx = round(linspace(8, nr-7, n_layers));
range_idx = reshape(layer_idx(randi(n_layers, n_scatter, 1)),[],1);
angle_idx = reshape(randi(na, n_scatter, 1),[],1);
amplitude_db = reshape(-dynamic_range_db * rand(n_scatter,1),[],1);
amplitude = 10.^(amplitude_db/20);
phase = 2*pi*rand(n_scatter,1);
coefficient = reshape(amplitude .* exp(1j*phase),[],1);

scene = complex(zeros(nr,na));
for k = 1:n_scatter
    scene(range_idx(k),angle_idx(k)) = ...
        scene(range_idx(k),angle_idx(k)) + coefficient(k);
end

range_m = reshape(cfg.range_axis_m(range_idx),[],1);
angle_deg = reshape(cfg.angle_grid_deg(angle_idx),[],1);
catalog = table(range_idx, angle_idx, range_m, angle_deg, amplitude_db, coefficient, ...
    'VariableNames', {'RangeIndex','AngleIndex','RangeM','AngleDeg', ...
    'AmplitudeDb','Coefficient'});
end
