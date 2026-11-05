function build_legacy_mex(backend)
% Build optional CUDA source against your supported MATLAB compiler/toolkit.
% This writes only to the local git-ignored build_mex directory.
if nargin < 1
    backend = "all";
end
backend = string(backend);
assert(isscalar(backend) && any(backend == ["echo","raytrace","all"]), ...
    'backend must be echo, raytrace, or all.');
root = legacy_sonarsim_setup();
assert(ispc || (isunix && ~ismac), 'This CUDA build recipe requires Windows or Linux with NVIDIA CUDA.');
out = fullfile(root, 'build_mex');
if ~isfolder(out)
    mkdir(out);
end
if any(backend == ["echo","all"])
    src = fullfile(root, 'sim_rx_from_scatterers_perTX_cuda_mex.cu');
    if ispc
        candidates = {fullfile(matlabroot, 'sys', 'cuda', 'win64', 'cuda', 'lib', 'cufft.lib')};
        toolkit = getenv('CUDA_PATH');
        if ~isempty(toolkit)
            candidates{end+1} = fullfile(toolkit, 'lib', 'x64', 'cufft.lib');
        end
        found = find(cellfun(@isfile, candidates), 1);
        assert(~isempty(found), 'cuFFT import library not found. Configure a supported MATLAB CUDA toolkit.');
        mexcuda('-R2018a', '-outdir', out, src, candidates{found});
    else
        mexcuda('-R2018a', '-outdir', out, src, '-lcufft');
    end
end
if any(backend == ["raytrace","all"])
    src = fullfile(root, 'Lib', 'stl_raytrace_arrivals_cuda3_5_1.cu');
    mexcuda('-R2018a', '-outdir', out, src);
end
addpath(out, '-begin');
end
