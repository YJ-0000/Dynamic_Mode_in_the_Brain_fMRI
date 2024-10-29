clear; clc;
current_path = pwd;

load('results/Schizo_timeseires_cortical_subcortical.mat');

data_path = 'G:\COBRE_preprocessed\4197885';

%%
% Define the frequency range for the band-pass filter
lowCutoff = 0.01;  % Hz
highCutoff = 0.1;  % Hz

%%
schizo_timeseries_extracted_denoised_filtered = cell(size(schizo_timeseries_extracted));
schizo_TR_info = cell(size(schizo_timeseries_extracted));
for nsub = 1:length(schizo_timeseries_extracted)
    sub_timeseries_extracted = schizo_timeseries_extracted{nsub};
    sub_name = schizo_sub_list{nsub};
    fprintf(['sub >> ',sub_name,'.\n']);

    tr = 2;

    Fs = 1/tr;

    fmri_data = sub_timeseries_extracted{1}';

    len_time = size(fmri_data,1);

    move_regressor = sub_timeseries_extracted{2}';
    move_regressor = [move_regressor, [zeros(1,6);(move_regressor(2:end,:)-move_regressor(1:end-1,:))]];
    move_regressor = [move_regressor, move_regressor.^2];

    noise_mat = [ones(len_time,1), (1:len_time)',((1:len_time).^2)',move_regressor, sub_timeseries_extracted{3}', sub_timeseries_extracted{4}'];

    fmri_data2 = (fmri_data - noise_mat*pinv(noise_mat)*(fmri_data));
    fmri_data2 = bandpass(fmri_data2, [lowCutoff,highCutoff], Fs)';

    fprintf('Filtered and aCompCor.\n');
    schizo_timeseries_extracted_denoised_filtered{nsub} = fmri_data2;
end

%%
save results/Schizo_timeseires_cortical_subcortical_denoised_filtered schizo_timeseries_extracted_denoised_filtered schizo_sub_list