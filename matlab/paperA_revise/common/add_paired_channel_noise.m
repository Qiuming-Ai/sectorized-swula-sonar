function [raw_noisy, diagnostics] = add_paired_channel_noise(raw_clean, snr_db, standardized_noise)
%ADD_PAIRED_CHANNEL_NOISE Add a prescribed standardized complex realization.
%   SNR is mean noise-free raw-channel signal power divided by mean added
%   raw-channel noise power.  The same standardized realization can be
%   reused across methods and then scaled to an equal per-method input SNR.

if ~isequal(size(raw_clean),size(standardized_noise))
    error('standardized_noise must have the same size as raw_clean.');
end
unit_noise = standardized_noise/sqrt(mean(abs(standardized_noise(:)).^2)+eps);
signal_power = mean(abs(raw_clean(:)).^2);
if isfinite(snr_db) && signal_power>0
    requested_noise_power = signal_power/10^(snr_db/10);
    added_noise = sqrt(requested_noise_power)*unit_noise;
else
    requested_noise_power = 0;
    added_noise = zeros(size(raw_clean),'like',raw_clean);
end
raw_noisy = raw_clean+added_noise;
actual_noise_power = mean(abs(added_noise(:)).^2);

diagnostics.signal_power = signal_power;
diagnostics.noise_power = actual_noise_power;
diagnostics.requested_snr_db = snr_db;
diagnostics.actual_snr_db = 10*log10(signal_power/(actual_noise_power+eps)+eps);
diagnostics.definition = "mean noise-free raw-channel signal power / mean added raw-channel noise power";
end
