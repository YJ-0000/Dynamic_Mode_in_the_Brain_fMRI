clear; clc;
task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};

for n_task = 1:length(task_names)
    current_task = task_names{n_task};
    
    load(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_meta.mat']);
    load(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted.mat']);
    load(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_task_input']);

    %%
    % Define the sampling frequency
    Fs = 1/0.72;  % Hz
    t_sample = 0.72;

    % Define the frequency range for the band-pass filter
    lowCutoff = 0.01;  % Hz
    highCutoff = 0.1;  % Hz

    %%
    time_series_preproc_filtered = cell(size(time_series_preproc));
    X = [ones(len_time,1), (1:len_time)',((1:len_time).^2)'];
    % X = [(1:len_time)',((1:len_time).^2)'];
    parfor nsub = 1:size(time_series_preproc,1)
        for nses = 1:2
            disp([nsub,nses]);
            data_mat = time_series_preproc{nsub,nses};
            task_info = task_input{nsub,nses};
            task_mat = [];
            M = length(task_info);
            for m = 1:M
                task_info1 = task_info{m};
                if ~isempty(task_info1)
                    onset_1= task_info1(:,1).Variables';
                    duration_1 = task_info1(:,2).Variables';
                    convolved_regressor1 = createTaskRegressor(t_sample, size(data_mat,2)*t_sample, onset_1, duration_1, 30);
                    task_mat = [task_mat,convolved_regressor1'];
                end
            end
            noise_mat = [noise_regressors{nsub,nses},X,task_mat];
            if ~isempty(data_mat)
                % apply residual matrix
                data_mat2 = (data_mat' - noise_mat*pinv(noise_mat)*(data_mat'));
                % bandpass filter
                data_mat2 = highpass(data_mat2, lowCutoff, Fs)';
                time_series_preproc_filtered{nsub,nses} = data_mat2;
            end
        end
    end

    %%
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_filtered_CompCor_meta'], 'sub_ids', 'folder_preproc_list', 'N', 'len_time')
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_filtered_CompCor.mat'], 'time_series_preproc_filtered', '-v7.3')
end