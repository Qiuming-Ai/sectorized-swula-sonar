function report = legacy_sonarsim_smoke(timing_dir)
% Read-only smoke test. Pass a directory with the exact timing inputs to
% exercise one CPU ping per method. No images or timing files are written.
if nargin < 1
    timing_dir = [];
end
root = legacy_sonarsim_setup();
assert(strcmpi(which('das_plane_wave_fd'), fullfile(root,'Lib','das_plane_wave_fd.m')), ...
    'Unexpected MATLAB path shadowing of the local FD beamformer.');
assert(nargout(@das_plane_wave_fd) == 2, 'Unexpected FD beamformer interface.');
report.recorded_mean_ms = verify_recorded_timing();
for file = ["PaperA_exp5_0.mat","PaperA_exp5_1.mat","PaperA_exp5_0_pc.mat","PaperA_exp7_0.mat"]
    info = whos('-file', legacy_figure_asset(char(file)));
    assert(~isempty(info), 'Legacy figure asset has no variables.');
end
rx = complex(zeros(32,2));
rx(4,:) = 1;
cfg = struct('t0',.04,'fs',125000,'c',1500,'rx_xyz',[0 0 0; .002 0 0], ...
    'angles',[-.1 0 .1],'f0',450000,'fd',struct('nfft',64,'sign',1));
[beam, ranges] = das_plane_wave_fd(cfg,rx);
assert(size(beam,1) == numel(ranges) && size(beam,2) == 3 && all(isfinite(beam(:))), ...
    'Legacy FD-DAS unit pulse did not produce a finite 3-angle output.');
report.fd_das_shape = size(beam);
report.one_ping_each = false;
if ~isempty(timing_dir)
    paths = legacy_desktop_timing_inputs(timing_dir);
    first_ping_s = [BaselineProcessFD(paths.dense,1), ...
        CDMProcessFD(paths.cdm,1), FDMProcessFD(paths.fdm,1,false)];
    assert(numel(first_ping_s) == 3 && all(isfinite(first_ping_s)) && all(first_ping_s > 0), ...
        'A one-ping timing smoke failed.');
    report.one_ping_each = true;
end
fprintf('Legacy SonarSim smoke passed. Optional exact-input one-ping check: %d.\n', report.one_ping_each);
end
