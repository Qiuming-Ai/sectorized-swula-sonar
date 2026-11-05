function image_complex = reconstruct_channel_image(cfg, raw_data, model)
%RECONSTRUCT_CHANNEL_IMAGE Matched filtering followed by common DAS imaging.

n_wave_sample = size(model.waveforms,1);
n_range = cfg.n_range_samples;
n_angle = numel(cfg.angle_grid_deg);
image_complex = complex(zeros(n_range,n_angle));
n_full = size(raw_data,1)+n_wave_sample-1;

for k = 1:model.n_waveforms
    matched = fft_filter_channels(raw_data,conj(flipud(model.waveforms(:,k))),n_full);
    range_data = matched(n_wave_sample:n_wave_sample+n_range-1,:);
    scan = model.output_sector{k};
    image_complex(:,scan) = image_complex(:,scan) + ...
        range_data*model.rx_weights{k}(:,scan);
end
end
