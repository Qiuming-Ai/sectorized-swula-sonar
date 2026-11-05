function export_review_imaging_data(sonarsim_dir)
%EXPORT_REVIEW_IMAGING_DATA Re-run the Paper-A imaging scripts and export data.
%
% This driver leaves the historical SonarSim scripts unchanged.  It reads each
% source script, applies the controlled revision settings in memory, executes
% it, and saves one MATLAB v7 file per method/SNR for the Python plotting stage.
% With no argument, use the repository's sibling legacy_sonarsim directory.
% An explicit argument must point to a trusted compatible SonarSim source set.
%
% Controlled revision settings:
%   * 100-degree FOV (-50 to 50 deg)
%   * six steering angles: +/-8, +/-24.25, +/-41.25 deg
%   * common 450-kHz wavelength for the physical aperture definitions
%   * one measured-AWGN realization at each declared input SNR
%   * deterministic random seeds
%   * a common continuous Hamming receive-aperture taper
%   * targets outside the revised FOV are excluded

this_dir = fileparts(mfilename('fullpath'));
workspace_root = fileparts(fileparts(this_dir));
if nargin < 1
    sonarsim_dir = fullfile(workspace_root, 'legacy_sonarsim');
    if ~isfolder(sonarsim_dir)
        sonarsim_dir = fullfile(workspace_root, 'SonarSim');
    end
end
assert(isfolder(sonarsim_dir), 'A8:MissingLegacySources', ...
    'A8 requires legacy_sonarsim or a trusted SonarSim directory: %s', sonarsim_dir);
[ok, source_info] = fileattrib(sonarsim_dir);
assert(ok, 'A8:UnreadableLegacySources', 'Cannot access %s.', sonarsim_dir);
sonarsim_dir = source_info.Name;
data_dir = fullfile(this_dir, 'data');

main_snrs = -15;
low_snrs = [-10, -20, -30];

jobs = [ ...
    make_job('main', 'dense',   'Dense ULA',          'PaperA_exp3_baseline_noise_doppler_baseband.m', 'beam',  'ranges_out', main_snrs, 1101), ...
    make_job('main', 'cdm',     'CDM SW-ULA',         'PaperA_exp3_cdm_noise_doppler_baseband.m',      'beams', 'ranges_out', main_snrs, 1201), ...
    make_job('main', 'fdm',     'FDM SW-ULA',         'PaperA_exp3_fdm_noise_doppler_baseband.m',      'beams', 'ranges_out', main_snrs, 1301), ...
    make_job('main', 'random',  'Random sparse',      'PaperA_exp4_random_noise_doppler_baseband.m',   'beam',  'ranges_out', main_snrs, 1401), ...
    make_job('main', 'nested',  'Nested sparse',      'PaperA_exp4_Nested_noise_doppler_baseband.m',   'beam',  'ranges_out', main_snrs, 1501), ...
    make_job('main', 'coprime', 'Coprime sparse',     'PaperA_exp4_Coprime_noise_doppler_baseband.m',  'beam',  'ranges_out', main_snrs, 1601), ...
    make_job('low_snr', 'dense', 'Dense ULA',         'PaperA_exp3_baseline_noise_doppler_baseband.m', 'beam',  'ranges_out', low_snrs, 2101), ...
    make_job('low_snr', 'cdm',   'CDM SW-ULA',        'PaperA_exp3_cdm_noise_doppler_baseband.m',      'beams', 'ranges_out', low_snrs, 2201), ...
    make_job('low_snr', 'fdm',   'FDM SW-ULA',        'PaperA_exp3_fdm_noise_doppler_baseband.m',      'beams', 'ranges_out', low_snrs, 2301) ...
    ];

source_names = unique({jobs.source_file});
for k = 1:numel(source_names)
    source_path = fullfile(sonarsim_dir, source_names{k});
    assert(isfile(source_path), 'A8:MissingLegacyEntry', ...
        'Missing historical A8 source script: %s', source_path);
end
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
old_path = path;
cleanup_path = onCleanup(@() path(old_path));
cd(sonarsim_dir);
addpath(genpath(sonarsim_dir));

old_visibility = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cleanup_fig = onCleanup(@() set(groot, 'defaultFigureVisible', old_visibility));

manifest_rows = cell(0, 9);
targets_written = false;

for job_idx = 1:numel(jobs)
    job = jobs(job_idx);
    fprintf('\n[%d/%d] %s / %s from %s\n', ...
        job_idx, numel(jobs), job.suite, job.display_name, job.source_file);

    expected_files = cell(1, numel(job.snr_db));
    for expected_idx = 1:numel(job.snr_db)
        expected_files{expected_idx} = fullfile(data_dir, sprintf( ...
            '%s_%s_snr_%s.mat', job.suite, job.method_key, ...
            snr_to_token(job.snr_db(expected_idx))));
    end
    if all(cellfun(@(path) exist(path, 'file') == 2, expected_files))
        fprintf('  skipped: all expected exports already exist\n');
        continue;
    end

    result = execute_patched_source(sonarsim_dir, job);

    if ~targets_written
        target_table = table(result.target_angle_deg(:), result.target_range_m(:), ...
            result.target_amplitude(:), ...
            'VariableNames', {'angle_deg', 'range_m', 'amplitude'});
        writetable(target_table, fullfile(data_dir, 'target_coordinates.csv'));
        targets_written = true;
    end

    for snr_idx = 1:numel(job.snr_db)
        image_linear = double(result.images{snr_idx, 1}); %#ok<NASGU>
        angle_deg = double(result.angles{snr_idx, 1}(:)); %#ok<NASGU>
        range_m = double(result.ranges{snr_idx, 1}(:)); %#ok<NASGU>
        target_angle_deg = double(result.target_angle_deg(:)); %#ok<NASGU>
        target_range_m = double(result.target_range_m(:)); %#ok<NASGU>
        target_amplitude = double(result.target_amplitude(:)); %#ok<NASGU>
        snr_db = double(job.snr_db(snr_idx)); %#ok<NASGU>
        random_seed = double(job.seed + snr_idx - 1); %#ok<NASGU>
        suite = job.suite; %#ok<NASGU>
        method_key = job.method_key; %#ok<NASGU>
        method_name = job.display_name; %#ok<NASGU>
        source_script = job.source_file; %#ok<NASGU>
        fov_deg = [-50, 50]; %#ok<NASGU>
        steering_deg = [-41.25, -24.25, -8, 8, 24.25, 41.25]; %#ok<NASGU>
        apodization = 'continuous Hamming taper over physical receive aperture'; %#ok<NASGU>
        noise_definition = 'single AWGN realization; SNR measured on each method input echo'; %#ok<NASGU>

        validate_export(image_linear, angle_deg, range_m, job);

        snr_token = snr_to_token(snr_db);
        out_name = sprintf('%s_%s_snr_%s.mat', suite, method_key, snr_token);
        out_path = fullfile(data_dir, out_name);
        save(out_path, 'image_linear', 'angle_deg', 'range_m', ...
            'target_angle_deg', 'target_range_m', 'target_amplitude', ...
            'snr_db', 'random_seed', 'suite', 'method_key', 'method_name', ...
            'source_script', 'fov_deg', 'steering_deg', 'apodization', ...
            'noise_definition', '-v7');

        manifest_rows(end+1, :) = {suite, method_key, method_name, snr_db, ...
            random_seed, job.source_file, size(image_linear, 1), ...
            size(image_linear, 2), out_name}; %#ok<AGROW>
        fprintf('  exported %s [%d x %d]\n', out_name, ...
            size(image_linear, 1), size(image_linear, 2));
    end
    close all force;
end

all_exports = dir(fullfile(data_dir, '*.mat'));
manifest_rows = cell(0, 9);
for export_idx = 1:numel(all_exports)
    export_path = fullfile(all_exports(export_idx).folder, all_exports(export_idx).name);
    meta = load(export_path, 'suite', 'method_key', 'method_name', 'snr_db', ...
        'random_seed', 'source_script');
    export_vars = whos('-file', export_path);
    image_info = export_vars(strcmp({export_vars.name}, 'image_linear'));
    manifest_rows(end+1, :) = {meta.suite, meta.method_key, meta.method_name, ...
        meta.snr_db, meta.random_seed, meta.source_script, ...
        image_info.size(1), image_info.size(2), all_exports(export_idx).name}; %#ok<AGROW>
end
manifest = cell2table(manifest_rows, 'VariableNames', ...
    {'suite', 'method_key', 'method_name', 'snr_db', 'random_seed', ...
     'source_script', 'n_range', 'n_angle', 'data_file'});
writetable(manifest, fullfile(data_dir, 'data_manifest.csv'));

fprintf('\nData export complete: %s\n', data_dir);
end


function job = make_job(suite, method_key, display_name, source_file, ...
    image_expression, range_expression, snr_db, seed)
job = struct( ...
    'suite', suite, ...
    'method_key', method_key, ...
    'display_name', display_name, ...
    'source_file', source_file, ...
    'image_expression', image_expression, ...
    'range_expression', range_expression, ...
    'snr_db', snr_db, ...
    'seed', seed);
end


function result = execute_patched_source(sonarsim_dir, job)
source_path = fullfile(sonarsim_dir, job.source_file);
assert(exist(source_path, 'file') == 2, 'Missing source script: %s', source_path);

code = fileread(source_path);
code = normalize_newlines(code);

% Do not allow a historical script to clear this function workspace.
code = regexprep(code, '(?m)^\s*clear\s+all\s*;?\s*close\s+all\s*;?\s*$', '');
code = regexprep(code, '(?m)^\s*close\s+all\s*;?\s*$', '');
code = regexprep(code, '(?m)^\s*figure\s*;?\s*$', '');

% Deterministic random geometry and noise.
code = sprintf('rng(%d, ''twister'');\n%s', job.seed, code);

% Use the requested, explicitly reported SNR values.
snr_literal = mat2str(job.snr_db);
code = regexprep(code, ...
    '(?m)^(?!\s*%)\s*snr_level\s*=\s*\[[^\]]*\]\s*;', ...
    sprintf('snr_level = %s;', snr_literal), 'once');

% Revised FOV and exact six-beam steering used by A1--A6.
code = strrep(code, ...
    'angles_deg = -[-47.5,-26,-8.5,8.5,26,47.5]'';', ...
    'angles_deg = -[-41.25,-24.25,-8,8,24.25,41.25]'';');
code = strrep(code, '(-60:0.1:60) * pi/180', '(-50:0.1:50) * pi/180');
code = strrep(code, '[-60:0.1:60]/180*pi', 'review_scan_deg/180*pi');

sector_code = [ ...
    'review_scan_deg = -50:0.1:50;' newline ...
    'review_edges_deg = [-50, -32.75, -16.125, 0, 16.125, 32.75, 50];' newline ...
    'angles_div = cell(1, 6);' newline ...
    'for review_sector_idx = 1:6' newline ...
    '    if review_sector_idx < 6' newline ...
    '        review_mask = review_scan_deg >= review_edges_deg(review_sector_idx) & review_scan_deg < review_edges_deg(review_sector_idx + 1);' newline ...
    '    else' newline ...
    '        review_mask = review_scan_deg >= review_edges_deg(review_sector_idx) & review_scan_deg <= review_edges_deg(review_sector_idx + 1);' newline ...
    '    end' newline ...
    '    angles_div{review_sector_idx} = review_scan_deg(review_mask);' newline ...
    'end'];
code = regexprep(code, ...
    '(?ms)^\s*angles_div\{1\}\s*=.*?^\s*angles_div\{6\}\s*=.*?;\s*$', ...
    sector_code, 'once');

% Ensure review_scan_deg exists for scripts without sector partitions.
code = regexprep(code, '(?m)^(\s*for\s+k\s*=\s*1:length\(snr_level\))', ...
    ['review_scan_deg = -50:0.1:50;' newline ...
     'review_images = cell(length(snr_level), length(velocity));' newline ...
     'review_ranges = cell(length(snr_level), length(velocity));' newline ...
     'review_angles = cell(length(snr_level), length(velocity));' newline ...
     '$1'], 'once');

% A single wavelength at the stated 450-kHz reference frequency gives a
% fixed receive aperture across Dense/CDM/FDM.
code = regexprep(code, '(?m)^\s*lambda\s*=\s*c0/f(?:min|max)\s*;', ...
    'lambda = c0/f0;');
code = strrep(code, ...
    'linspace(-probe.N*probe.pitch/2,probe.N*probe.pitch/2,probe.N)', ...
    'linspace(-(probe.N-1)*probe.pitch/2,(probe.N-1)*probe.pitch/2,probe.N)');
code = strrep(code, ...
    'linspace(-probe.N0*probe.pitch/2,probe.N0*probe.pitch/2,probe.N0)', ...
    'linspace(-(probe.N0-1)*probe.pitch/2,(probe.N0-1)*probe.pitch/2,probe.N0)');

% Remove the +/-55-degree targets that lie outside the revised FOV.
target_filter = [ ...
    'review_target_angle_deg = atan2d(point_position(:,1), point_position(:,3));' newline ...
    'review_target_keep = abs(review_target_angle_deg) <= 50 + 1e-9;' newline ...
    'point_position = point_position(review_target_keep, :);' newline ...
    '% Set point amplitudes'];
code = regexprep(code, '(?m)^% Set point amplitudes', target_filter, 'once');

% Use the same continuous Hamming taper as a function of physical position.
apodization_code = [ ...
    'review_x = rx_xyz(:,1);' newline ...
    'review_x_norm = (review_x - min(review_x)) ./ max(max(review_x) - min(review_x), eps);' newline ...
    'win = 0.54 - 0.46*cos(2*pi*review_x_norm);'];
code = regexprep(code, '(?m)^(?!\s*%)\s*win\s*=[^\r\n]*$', ...
    apodization_code, 'once');

% Historical CDM/FDM and non-uniform scripts sum four independent noisy
% copies, which changes the declared SNR.  Replace that block by one measured
% AWGN realization so every exported panel has the same SNR definition.
code = regexprep(code, ...
    ['(?ms)^\s*v_calc_1\s*=\s*awgn\(v_calc,snr_level\(k\),''measured''\);\s*' ...
     '^\s*v_calc_2\s*=\s*awgn\(v_calc,snr_level\(k\),''measured''\);\s*' ...
     '^\s*v_calc_3\s*=\s*awgn\(v_calc,snr_level\(k\),''measured''\);\s*' ...
     '^\s*v_calc_4\s*=\s*awgn\(v_calc,snr_level\(k\),''measured''\);\s*' ...
     '^\s*v_calc\s*=\s*v_calc_1\+v_calc_2\+v_calc_3\+v_calc_4\s*;'], ...
    '        rng(review_seed_base + k - 1, ''twister'');\n        v_calc = awgn(v_calc, snr_level(k), ''measured'');');

% Apply a deterministic per-SNR seed to scripts that already use one AWGN call.
code = regexprep(code, ...
    '(?m)^(\s*)v_calc\s*=\s*awgn\(v_calc,snr_level\(k\),''measured''\);', ...
    ['$1rng(review_seed_base + k - 1, ''twister'');' newline ...
     '$1v_calc = awgn(v_calc, snr_level(k), ''measured'');']);
code = sprintf('review_seed_base = %d;\n%s', job.seed, code);

% Capture numerical arrays instead of relying on MATLAB figure objects.
capture_code = sprintf([ ...
    '        review_images{k,j} = abs(%s);' newline ...
    '        review_ranges{k,j} = %s;' newline ...
    '        review_angles{k,j} = cfg.angles(:);'], ...
    job.image_expression, job.range_expression);
code = regexprep(code, '(?m)^\s*sector_plot\([^\r\n]*$', capture_code);

% Remove plot-only commands left after capture.  All publication rendering is
% performed by Python from the exported arrays.
code = regexprep(code, '(?m)^\s*ax\s*=\s*subplot\([^\r\n]*$', '');
code = regexprep(code, '(?m)^\s*hold\s+on\s*;?\s*$', '');
code = regexprep(code, '(?m)^\s*scatter3\([^\r\n]*$', '');
code = regexprep(code, '(?m)^\s*(xlim|ylim|view|title)\([^\r\n]*$', '');

% Keep the exact controlled script for provenance and troubleshooting.
debug_dir = fullfile(fileparts(mfilename('fullpath')), 'debug_patched');
if ~exist(debug_dir, 'dir')
    mkdir(debug_dir);
end
debug_name = sprintf('%s_%s_%s', job.suite, job.method_key, job.source_file);
debug_path = fullfile(debug_dir, debug_name);
debug_fid = fopen(debug_path, 'w');
assert(debug_fid >= 0, 'Cannot create debug script: %s', debug_path);
debug_cleanup = onCleanup(@() fclose(debug_fid));
fwrite(debug_fid, code, 'char');
clear debug_cleanup;

% Run the controlled in-memory script.
eval(code);

assert(exist('review_images', 'var') == 1, 'No image data captured from %s', job.source_file);
assert(all(~cellfun(@isempty, review_images(:,1))), ...
    'One or more SNR images were not captured from %s', job.source_file);

target_angle_deg = atan2d(point_position(:,1), point_position(:,3));
target_range_m = hypot(point_position(:,1), point_position(:,3));

result = struct;
result.images = review_images;
result.ranges = review_ranges;
result.angles = cellfun(@(x) rad2deg(x(:)), review_angles, 'UniformOutput', false);
result.target_angle_deg = target_angle_deg;
result.target_range_m = target_range_m;
result.target_amplitude = point_amplitudes;
end


function validate_export(image_linear, angle_deg, range_m, job)
assert(isnumeric(image_linear) && ismatrix(image_linear), ...
    '%s did not produce a numeric 2-D image.', job.source_file);
assert(all(isfinite(image_linear(:))), ...
    '%s produced non-finite image values.', job.source_file);
assert(numel(angle_deg) == size(image_linear, 2), ...
    '%s angle dimension mismatch: %d labels versus %d columns.', ...
    job.source_file, numel(angle_deg), size(image_linear, 2));
assert(numel(range_m) == size(image_linear, 1), ...
    '%s range dimension mismatch: %d labels versus %d rows.', ...
    job.source_file, numel(range_m), size(image_linear, 1));
assert(min(angle_deg) >= -50.05 && max(angle_deg) <= 50.05, ...
    '%s exported angles outside the revised FOV.', job.source_file);
end


function token = snr_to_token(snr_db)
if snr_db < 0
    token = sprintf('m%g', abs(snr_db));
else
    token = sprintf('p%g', snr_db);
end
token = strrep(token, '.', 'p');
end


function text = normalize_newlines(text)
text = strrep(text, sprintf('\r\n'), sprintf('\n'));
text = strrep(text, sprintf('\r'), sprintf('\n'));
end
