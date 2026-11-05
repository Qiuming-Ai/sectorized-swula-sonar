function timing_table = CollectTimingData(input_dir, output_csv, show_plots)
% New replay of the historical Dense -> FDM -> CDM desktop workloads.
% Optional input directory defaults to data/desktop-timing/overlay.
% No file is written unless an explicit, unused output CSV is provided.
legacy_sonarsim_setup();
if nargin < 1, input_dir = []; end
if nargin < 2, output_csv = []; end
if nargin < 3, show_plots = false; end
paths = legacy_desktop_timing_inputs(input_dir);
files = struct2cell(paths);
assert(all(cellfun(@isfile, files)), 'All three exact synthetic timing HDF5 inputs are required.');
if ~isempty(output_csv)
    assert(~isfile(output_csv), 'Refusing to overwrite timing CSV: %s', output_csv);
    parent = fileparts(output_csv);
    assert(isempty(parent) || isfolder(parent), 'Create the output directory before collecting timing.');
end
baseline_timing = BaselineProcessFD(paths.dense);
fdm_timing = FDMProcessFD(paths.fdm, Inf, show_plots);
cdm_timing = CDMProcessFD(paths.cdm);
assert(all([numel(baseline_timing), numel(fdm_timing), numel(cdm_timing)] == 100), ...
    'Expected the complete 100-ping historical fixture for each method.');
timing_table = table(baseline_timing, fdm_timing, cdm_timing, ...
    'VariableNames', {'Baseline_Time','FDM_Time','CDM_Time'});
if ~isempty(output_csv)
    writetable(timing_table, output_csv);
end
fprintf('New replay, mean ms (Dense / CDM / FDM): %.5f / %.5f / %.5f\n', ...
    1000 * mean(baseline_timing), 1000 * mean(cdm_timing), 1000 * mean(fdm_timing));
end
