function binary = build_shipwreck_mex(out_dir)
%BUILD_SHIPWRECK_MEX Compile the released CUDA source with local MATLAB.
% Requires mexcuda, a supported C++ toolchain and CUDA setup for this MATLAB.
here = fileparts(mfilename('fullpath'));
package_root = fileparts(fileparts(here));
if nargin<1 || isempty(out_dir)
    out_dir = fullfile(package_root,'outputs','matlab','shipwreck_cuda');
end
if exist('mexcuda','file')==0
    error('Shipwreck:MissingMexcuda','MATLAB mexcuda is required; see shipwreck_cuda/README.md.');
end
if ~isfolder(out_dir), mkdir(out_dir); end
source = fullfile(here,'stl_raytrace_arrivals_cuda3_5_1.cu');
mexcuda('-R2018a','-outdir',out_dir,source);
binary = fullfile(out_dir,['stl_raytrace_arrivals_cuda3_5_1.' mexext]);
assert(isfile(binary),'mexcuda did not produce the expected binary: %s',binary);
fprintf('Built local MEX: %s\n',binary);
end
