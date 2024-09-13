function convolved_regressor = createTaskRegressor(TR, total_time, onsets, durations, upsample_factor)
    % Function to create a task regressor convoluted with an HRF
    % Inputs:
    %   TR - Time of repetition in seconds
    %   total_time - Total duration of the experiment in seconds
    %   onsets - Array of onset times of the task in seconds
    %   durations - Array of durations of each task event in seconds
    %   upsample_factor - Factor to upsample the original time series
    
    % Create a time series vector with 1s at the onset of each task event
    t = 0:(TR/upsample_factor):total_time;
    upsampled_time_series = zeros(1, numel(t));
    for i = 1:length(onsets)
        target_idx = (onsets(i) <= t) & ((onsets(i)+durations(i)) >= t);
        upsampled_time_series(target_idx) = 1;
    end

    % Upsample the time series
%     time_series = repelem(upsampled_time_series, 1/upsample_factor);

    % Convolve the upsampled time series with the canonical HRF
    hrf = spm_hrf(TR / upsample_factor);  % Adjust HRF calculation for upsampled TR
    convolved_regressor_full = conv(upsampled_time_series, hrf);
    
    % Downsample the convolved regressor back to the original TR
    convolved_regressor = convolved_regressor_full(ceil(upsample_factor/2):upsample_factor:end);
    convolved_regressor = convolved_regressor(1:ceil(total_time / TR));  % Ensure it matches the experiment's duration

    % Return the convolved regressor
%     t_down = t(ceil(upsample_factor/2):upsample_factor:end); t_down = t_down(1:ceil(total_time / TR)); 
%     figure; scatter(t_down,convolved_regressor); hold on; plot(t,upsampled_time_series);
end
