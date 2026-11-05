function path = legacy_model_path(filename)
% Provide the exact optional original mesh locally, subject to its rights.
root = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(root));
if strcmp(filename, 'ship.stl')
    extracted = fullfile(repo, 'data', 'shipwreck', 'geometry', filename);
    if isfile(extracted)
        path = extracted;
        return;
    end
end
path = fullfile(root, 'Models', filename);
assert(isfile(path), ['Legacy Exp8 requires local mesh %s. Models are not ' ...
    'distributed in this repository. See matlab/legacy_sonarsim/README.md.'], path);
end
