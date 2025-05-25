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
    
    load secure_data/path_info;
    HCP_preprocessed_path = [HCP_preprocessed_task_path, current_task];
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
    
    %% Identifying problematic subjects
    cd(current_path);
    load secure_data/path_info;
    behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
    for nrow = size(behav_data_table,1):-1:1
        if ~any(sub_ids==behav_data_table(nrow,'Subject').Variables)
            behav_data_table(nrow,:) = [];
        end
    end
    behav_data_table = sortrows(behav_data_table, 'Subject');

    does_have_MMSE = ~strcmp(behav_data_table.MMSE_Compl,'true');
    MMSE_thres = 26;
    is_cognitive_impaired = behav_data_table.MMSE_Score <= MMSE_thres;

    HCP_preproc_dir = folder_preproc_list;

    sub_ids_foler = zeros(size(HCP_preproc_dir));
    movement_thres = 0.15;
    is_excluded_due_movement = false(length(sub_ids),1);
    for n_fol = 1:length(HCP_preproc_dir)
        fol_prerpoc_name = HCP_preproc_dir(n_fol).name;
        aa = split(fol_prerpoc_name,'_');
        sub_idx = find(sub_ids == str2double(aa{1}));

        if isempty(sub_idx)
            continue;
        end

        fprintf('Checking %s whether mean relative movement exceeds %0.4f mm \n',aa{1},movement_thres)

        REST_1_or_2 = aa{4};
        rest_idx = 1;
        try
            movement_rel_RMS_LR_path = fullfile(HCP_preprocessed_path,HCP_preproc_dir(n_fol).name,aa{1},...
                'MNINonLinear','Results',['tfMRI_',REST_1_or_2,'_LR'],'Movement_RelativeRMS_mean.txt');
            movement_rel_RMS_LR = readmatrix(movement_rel_RMS_LR_path);
            movement_rel_RMS_RL_path = fullfile(HCP_preprocessed_path,HCP_preproc_dir(n_fol).name,aa{1},...
                'MNINonLinear','Results',['tfMRI_',REST_1_or_2,'_RL'],'Movement_RelativeRMS_mean.txt');
            movement_rel_RMS_RL = readmatrix(movement_rel_RMS_RL_path);
            is_excluded_due_movement(sub_idx,rest_idx) = ...
                any(movement_rel_RMS_LR > movement_thres) || any(movement_rel_RMS_RL > movement_thres);
        catch
            is_excluded_due_movement(sub_idx,rest_idx) = true;
        end
    end
    
    if strcmp(current_task,'WM')
        RL_processing_errors_subs = [168139;196952];
    elseif strcmp(current_task,'EMOTION')
        RL_processing_errors_subs = [168139;186545;320826;644044];
    elseif strcmp(current_task,'MOTOR')
        RL_processing_errors_subs = [144428;168139;192237;822244;870861;947668];
    elseif strcmp(current_task,'LANGUAGE')
        RL_processing_errors_subs = [168139]; %#ok<NBRAK>
    elseif strcmp(current_task,'GAMBLING')
        RL_processing_errors_subs = [168139;144428;822244;870861;947668];
    elseif strcmp(current_task,'SOCIAL')
        RL_processing_errors_subs = [168139;186545;748662;809252];
    elseif strcmp(current_task,'RELATIONAL')
        RL_processing_errors_subs = [168139;186545;223929;644044];
    
    end

    is_RL_processing_errors = false(length(sub_ids),1);
    for n_rl = 1:length(RL_processing_errors_subs)
        is_RL_processing_errors(sub_ids==RL_processing_errors_subs(n_rl)) = true;
    end

    %%
    save(['results/HCP_timeseries_tfMRI_',current_task,'_subject_exclude_info'], 'sub_ids', 'is_excluded_due_movement', 'does_have_MMSE', 'is_cognitive_impaired', 'MMSE_thres', 'movement_thres', 'is_RL_processing_errors');
    
end