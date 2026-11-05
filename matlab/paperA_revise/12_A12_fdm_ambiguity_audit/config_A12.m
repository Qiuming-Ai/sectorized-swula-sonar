function cfg = config_A12()
%CONFIG_A12 Full-delay ambiguity audit for the implemented waveform bank.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));
cfg = exp_default_config("paper");
cfg.methods = ["cdm","fdm_current"];
cfg.method_labels = ["CDM--ZC","packed FDM--LFM"];
cfg.velocities_mps = [-2,-1,-0.5,0,0.5,1,2];
cfg.processing_half_window_samples = 240;
cfg.full_delay_definition = "all finite-pulse lags from -(N-1) to +(N-1) samples";
cfg.processing_window_definition = "absolute lag no larger than 240 samples (0.6 ms at 400 kHz)";
cfg.zero_delay_definition = "lag equal to zero samples";
cfg.spectrum_fft_length = 2^nextpow2(16*round(cfg.pulse_duration_s*cfg.waveform_fs));
end
