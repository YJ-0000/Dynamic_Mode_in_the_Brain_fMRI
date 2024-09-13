clear;clc;
current_path = pwd;
%%

task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};
len_time_list = [405,176,284,316,253,274,232];

for n_task = 1:length(task_names)
    current_task = task_names{n_task};
    
    labels_ref = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
    labels_cortical = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
    labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
    label_cortical_data = round(labels_cortical.cdata);

    try
        label_subcortical_data = labels_subcortical.cdata;
        label_data = [label_cortical_data; label_subcortical_data(length(label_cortical_data)+1:end)];
    catch
        label_data = label_cortical_data;
    end

    HCP_preprocessed_path = ['G:\HCP_3T_', current_task];
    cd(HCP_preprocessed_path);
    folder_preproc_list = dir('*_preproc');
    folder_preproc_list = folder_preproc_list([folder_preproc_list.isdir]);

    %% sub IDs
    sub_ids_foler = zeros(size(folder_preproc_list));
    for n_fol = 1:length(folder_preproc_list)
        fol_preproc_name = folder_preproc_list(n_fol).name;
        aa = split(fol_preproc_name,'_');
        sub_ids_foler(n_fol) = str2double(aa{1});
    end

    sub_ids = unique(sub_ids_foler);
    num_rest_session = zeros(size(sub_ids));
    for n_sub = 1:length(sub_ids)
        num_rest_session(n_sub) = sum(sub_ids_foler==sub_ids(n_sub));
    end

    time_series_preproc = cell(length(sub_ids),2);
    noise_regressors = cell(length(sub_ids),2);
    %% ROI time-series
    len_time = len_time_list(n_task);

    label_idx_list = sort(unique(label_data));
    label_idx_list(label_idx_list==0) = [];
    N = numel(label_idx_list);

    for n_fol = 1:length(folder_preproc_list)
        cd(HCP_preprocessed_path);
        fol_preproc_name = folder_preproc_list(n_fol).name;
        cd(fol_preproc_name);
        aa = split(fol_preproc_name,'_'); 

        fol_preproc_name = folder_preproc_list(n_fol).name;
        bb = split(fol_preproc_name,'_'); 

        if ~strcmp(aa{1},bb{1}) || ~strcmp(aa{4},bb{4})
            error('Name mismatch');
        end

        cd(aa{1});
        cd('MNINonLinear\Results');

        sub_idx = find(sub_ids == str2double(aa{1}));

        LR_dir = dir('tfMRI_*_*');

        if length(LR_dir) < 2
            continue
        end

        cd(['tfMRI_',current_task,'_LR']);
        dtseries_file = dir('tfMRI_*_Atlas_MSMAll.dtseries.nii');
        if length(dtseries_file) < 1; continue; end
        cifti = cifti_read(dtseries_file(1).name);
    %     wm_regressor1 = readmatrix('rfMRI_REST1_LR_WM.txt');
    %     csf_regressor1 = readmatrix('rfMRI_REST1_LR_CSF.txt');
        try
            move_regressor1 = readmatrix('Movement_Regressors.txt');
        catch
            continue
        end

        if isempty(time_series_preproc{sub_idx,1})
            empty_col = 1;
        else
            empty_col = 3;
        end

        if size(cifti.cdata,2) ~= len_time
            continue
        end

        time_series_single1 = zeros(N,len_time);
        cifti_data = cifti.cdata(1:length(label_data),:);
        for n_region = 1:N
            time_series_single1(n_region,:) = mean(cifti_data(label_data==label_idx_list(n_region),:));
        end

        cd(HCP_preprocessed_path);
        cd(fol_preproc_name);
        cd(aa{1});
        cd('MNINonLinear\Results');


        cd(['tfMRI_',current_task,'_RL']);
        dtseries_file = dir('tfMRI_*_Atlas_MSMAll.dtseries.nii');
        if length(dtseries_file) < 1; continue; end
        cifti = cifti_read(dtseries_file(1).name);
    %     wm_regressor1 = readmatrix('rfMRI_REST1_LR_WM.txt');
    %     csf_regressor1 = readmatrix('rfMRI_REST1_LR_CSF.txt');
        try
            move_regressor2 = readmatrix('Movement_Regressors.txt');
        catch
            continue
        end

        if size(cifti.cdata,2) ~= len_time
            continue
        end

        time_series_single2 = zeros(N,len_time);
        cifti_data = cifti.cdata(1:length(label_data),:);
        for n_region = 1:N
            time_series_single2(n_region,:) = mean(cifti_data(label_data==label_idx_list(n_region),:));
        end

        time_series_preproc{sub_idx,empty_col} = time_series_single1;
        noise_regressors{sub_idx,empty_col} = [move_regressor1,move_regressor1.^2];
        empty_col = empty_col + 1;
        time_series_preproc{sub_idx,empty_col} = time_series_single2;
        noise_regressors{sub_idx,empty_col} = [move_regressor2,move_regressor2.^2];

        disp([n_fol,sub_ids(sub_idx),empty_col-1,empty_col]);
    end


    %%
    cd(current_path);
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted_meta'], 'sub_ids', 'folder_preproc_list', 'N', 'len_time', 'noise_regressors', 'label_idx_list');
    save(['results/HCP_timeseries_tfMRI_',current_task,'_cortical_subcortical_extracted.mat'], 'time_series_preproc', '-v7.3')
end