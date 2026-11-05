function root = legacy_sonarsim_setup()
% Add only this archive's helpers, not other workspace projects/toolboxes.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'Lib'), '-begin');
build_dir = fullfile(root, 'build_mex');
if isfolder(build_dir)
    addpath(build_dir, '-begin');
end
end
