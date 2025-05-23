clear; clc;
load('results/HCP_timeseries_cortical_subcortical_extracted_meta.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted.mat');

%%
% Define the sampling frequency
Fs = 1/0.72;  % Hz

% Define the frequency range for the band-pass filter
lowCutoff = 0.01;  % Hz
highCutoff = 0.1;  % Hz


x = time_series_denoised{1,1};
% Apply the filter
filteredData = bandpass(x(1,:),[lowCutoff, highCutoff],Fs);

%%
time_series_denoised_filtered = cell(size(time_series_denoised));
for nsub = 1:size(time_series_denoised,1)
    for nses = 1:4
        disp([nsub,nses]);
        data_mat = time_series_denoised{nsub,nses};
        if ~isempty(data_mat)
            data_mat = bandpass(data_mat',[lowCutoff, highCutoff],Fs)';
            time_series_denoised_filtered{nsub,nses} = data_mat;
        end
    end
end

%%
save results/HCP_timeseries_cortical_subcortical_extracted_filtered_meta sub_ids folder_denoised_list N lowCutoff highCutoff
save('results/HCP_timeseries_cortical_subcortical_extracted_filtered.mat', 'time_series_denoised_filtered', '-v7.3')