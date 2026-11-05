function paths = legacy_desktop_timing_inputs(input_dir)
% Exact historical synthetic input filenames; no substitutions.
if nargin < 1 || isempty(input_dir)
    root = fileparts(mfilename('fullpath'));
    repo = fileparts(fileparts(root));
    input_dir = fullfile(repo, 'data', 'desktop-timing', 'overlay');
end
paths = struct( ...
    'dense', fullfile(input_dir, 'BaselineData_20251121_150913.h5'), ...
    'cdm', fullfile(input_dir, 'CDMData_20251125_160029.h5'), ...
    'fdm', fullfile(input_dir, 'FDMData_20251126_112556.h5'));
end
