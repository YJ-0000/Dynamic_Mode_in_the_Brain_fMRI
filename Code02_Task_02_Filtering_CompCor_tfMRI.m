clear; clc;
task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};

for n_task = 1:length(task_names)
    current_task = task_names{n_task};
    
    load(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_meta.mat']);
    load(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted.mat']);

    %%
    % Define the sampling frequency
    Fs = 1/0.72;  % Hz

    % Define the frequency range for the band-pass filter
    lowCutoff = 0.01;  % Hz
    highCutoff = 0.1;  % Hz


    x = time_series_preproc{1,1};
    % Apply the filter
    filteredData = bandpass(x(1,:),[lowCutoff, highCutoff],Fs);

    %%
    time_series_preproc_filtered = cell(size(time_series_preproc));
    X = [ones(len_time,1), (1:len_time)',((1:len_time).^2)'];
    % X = [(1:len_time)',((1:len_time).^2)'];
    parfor nsub = 1:size(time_series_preproc,1)
        for nses = 1:2
            disp([nsub,nses]);
            data_mat = time_series_preproc{nsub,nses};
            noise_mat = [noise_regressors{nsub,nses},X];
            if ~isempty(data_mat)
                for n_roi = 1:N
                    % apply residual matrix
                    data_mat(n_roi,:) = data_mat(n_roi,:)' - noise_mat*(noise_mat\data_mat(n_roi,:)');
                    % bandpass filter
                    data_mat(n_roi,:) = ...
                        bandpass(data_mat(n_roi,:),[lowCutoff, highCutoff],Fs);
                end
                time_series_preproc_filtered{nsub,nses} = data_mat;
            end
        end
    end

    %%
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_filtered_CompCor_meta'], 'sub_ids', 'folder_preproc_list', 'N', 'len_time')
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_filtered_CompCor.mat'], 'time_series_preproc_filtered', '-v7.3')
end