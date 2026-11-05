function mean_ms = verify_recorded_timing()
% Reproduce the published desktop means without HDF5 inputs or new timing.
root = fileparts(mfilename('fullpath'));
T = readtable(fullfile(root, 'data', 'TimingData_20251126_114034.csv'));
assert(height(T) == 100 && isequal(T.Properties.VariableNames, ...
    {'Baseline_Time','FDM_Time','CDM_Time'}), 'Unexpected historical timing CSV schema.');
values = T{:,{'Baseline_Time','CDM_Time','FDM_Time'}};
assert(all(isfinite(values(:))) && all(values(:)>0), 'Invalid recorded timing.');
mean_ms = 1000 * mean(values,1);
assert(isequal(round(mean_ms,2), [265.61 60.11 31.27]), 'Reported means do not match the manuscript.');
fprintf('Verified 100 recorded calls/method: Dense %.6f, CDM %.6f, FDM %.6f ms.\n', mean_ms);
end
