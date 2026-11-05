function [image_complex, diagnostics] = simulate_reduced_image(cfg, scene, method, snr_db, seed)
%SIMULATE_REDUCED_IMAGE Reduced-order angle-range forward/adjoint model.
% This model preserves array responses, waveform cross-correlations and
% transmit-slot separation. Noise
% is applied at the reconstructed-output level; raw-channel claims must still
% be verified with the full SonarSim pipeline before manuscript submission.

if nargin < 4 || isempty(snr_db)
    snr_db = Inf;
end
if nargin < 5
    seed = cfg.default_seed;
end
model = prepare_reduced_model(cfg,method);
[image_complex,diagnostics] = apply_reduced_model(cfg,scene,model,snr_db,seed);
end
