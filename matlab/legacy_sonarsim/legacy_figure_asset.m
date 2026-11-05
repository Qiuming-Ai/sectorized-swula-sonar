function path = legacy_figure_asset(filename)
% The archived small MAT files are inputs, not locations for new exports.
root = fileparts(mfilename('fullpath'));
path = fullfile(root, 'data', 'figures', filename);
assert(isfile(path), 'Legacy figure asset is missing: %s', path);
end
