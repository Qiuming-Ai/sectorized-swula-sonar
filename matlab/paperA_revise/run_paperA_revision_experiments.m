function outputs = run_paperA_revision_experiments(mode,selection)
%RUN_PAPERA_REVISION_EXPERIMENTS Run numbered revision experiments in order.
% Examples:
%   run_paperA_revision_experiments("quick",1:3)
%   run_paperA_revision_experiments("paper",1:7)
if nargin < 1, mode = "quick"; end
if nargin < 2, selection = 1:7; end
root = fileparts(mfilename('fullpath'));
folders = ["01_A1_tx_beampattern","02_A3_cdm_scene_complexity", ...
    "03_A4_fdm_band_config","04_A6_snr_monte_carlo", ...
    "05_A2_two_way_psf","06_A5_doppler_sensitivity", ...
    "07_A7_baseline_fairness"];
functions = ["run_A1_tx_beampattern","run_A3_cdm_scene_complexity", ...
    "run_A4_fdm_band_config","run_A6_snr_monte_carlo", ...
    "run_A2_two_way_psf","run_A5_doppler_sensitivity", ...
    "run_A7_baseline_fairness"];
if any(selection<1 | selection>7 | fix(selection)~=selection)
    error('selection must contain integers from 1 through 7.');
end
addpath(fullfile(root,'common'));
for k=1:numel(folders), addpath(fullfile(root,folders(k))); end

outputs = struct();
for k=selection(:).'
    fprintf('\n===== Starting %02d/07: %s (%s mode) =====\n', ...
        k,folders(k),mode);
    timer=tic;
    outputs.(sprintf('Experiment%02d',k)) = feval(functions(k),mode);
    fprintf('===== Finished %02d/07 in %.1f s =====\n',k,toc(timer));
end
end
