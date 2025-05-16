clear;clc;
current_path = pwd;
%%
labels_ref = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
labels_cortical = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
label_cortical_data = round(labels_cortical.cdata);
labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');

try
    label_subcortical_data = round(labels_subcortical.cdata);
    label_data = [label_cortical_data; label_subcortical_data(length(label_cortical_data)+1:end)];
catch
    label_data = label_cortical_data;
end

load secure_data/path_info;
cd(HCP_denoised_path);
folder_denoised_list = dir('*_fixextended');
folder_denoised_list = folder_denoised_list([folder_denoised_list.isdir]);

%% sub IDs
sub_ids_foler = zeros(size(folder_denoised_list));
for n_fol = 1:length(folder_denoised_list)
    fol_denoised_name = folder_denoised_list(n_fol).name;
    aa = split(fol_denoised_name,'_');
    sub_ids_foler(n_fol) = str2double(aa{1});
end

sub_ids = unique(sub_ids_foler);
num_rest_session = zeros(size(sub_ids));
for n_sub = 1:length(sub_ids)
    num_rest_session(n_sub) = sum(sub_ids_foler==sub_ids(n_sub));
end

time_series_denoised = cell(length(sub_ids),4);
noise_regressors = cell(length(sub_ids),4);
%% ROI time-series
len_time = 1200;

label_idx_list = sort(unique(label_data));
label_idx_list(label_idx_list==0) = [];
N = numel(label_idx_list);

for n_fol = 1:length(folder_denoised_list)
    cd(HCP_denoised_path);
    fol_denoised_name = folder_denoised_list(n_fol).name;
    cd(fol_denoised_name);
    aa = split(fol_denoised_name,'_'); 
    
    cd(aa{1});
    cd('MNINonLinear\Results');
    
    sub_idx = find(sub_ids == str2double(aa{1}));
        
    LR_dir = dir('rfMRI_REST*');
    
    if length(LR_dir) < 2
        continue
    end
    
    try 
        cd('rfMRI_REST1_LR');
        cifti = cifti_read('rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii');
    catch
        cd('rfMRI_REST2_LR');
        cifti = cifti_read('rfMRI_REST2_LR_Atlas_hp2000_clean.dtseries.nii');
    end
    
    if isempty(time_series_denoised{sub_idx,1})
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
    
    cd(HCP_denoised_path);
    cd(fol_denoised_name);
    cd(aa{1});
    cd('MNINonLinear\Results');
    try
        cd('rfMRI_REST1_RL');
        cifti = cifti_read('rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii');
    catch
        cd('rfMRI_REST2_RL');
        cifti = cifti_read('rfMRI_REST2_RL_Atlas_hp2000_clean.dtseries.nii');
    end
    
    if size(cifti.cdata,2) ~= len_time
        continue
    end
    
    time_series_single2 = zeros(N,len_time);
    cifti_data = cifti.cdata(1:length(label_data),:);
    for n_region = 1:N
        time_series_single2(n_region,:) = mean(cifti_data(label_data==label_idx_list(n_region),:));
    end
    
    time_series_denoised{sub_idx,empty_col} = time_series_single1;
    empty_col = empty_col + 1;
    time_series_denoised{sub_idx,empty_col} = time_series_single2;
        
    disp([n_fol,sub_ids(sub_idx),empty_col-1,empty_col]);
end


%%
cd(current_path);
save results/HCP_timeseries_cortical_subcortical_extracted_meta sub_ids folder_denoised_list N len_time label_idx_list
save('results/HCP_timeseries_cortical_subcortical_extracted.mat', 'time_series_denoised', '-v7.3')

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

HCP_preproc_dir = dir(fullfile(HCP_preprocessed_rest_path,'*rfMRI*'));

sub_ids_foler = zeros(size(HCP_preproc_dir));
movement_thres = 0.15;
is_excluded_due_movement = false(length(sub_ids),2);
for n_fol = 1:length(HCP_preproc_dir)
    fol_prerpoc_name = HCP_preproc_dir(n_fol).name;
    aa = split(fol_prerpoc_name,'_');
    sub_idx = find(sub_ids == str2double(aa{1}));
    
    if isempty(sub_idx)
        continue;
    end
    
    fprintf('Checking %s whether mean relative movement exceeds %0.4f mm \n',aa{1},movement_thres)
    
    REST_1_or_2 = aa{4};
    rest_idx = str2double(REST_1_or_2(end));
    try
        movement_rel_RMS_LR_path = fullfile(HCP_preprocessed_rest_path,HCP_preproc_dir(n_fol).name,aa{1},...
            'MNINonLinear','Results',['rfMRI_',REST_1_or_2,'_LR'],'Movement_RelativeRMS_mean.txt');
        movement_rel_RMS_LR = readmatrix(movement_rel_RMS_LR_path);
        movement_rel_RMS_RL_path = fullfile(HCP_preprocessed_rest_path,HCP_preproc_dir(n_fol).name,aa{1},...
            'MNINonLinear','Results',['rfMRI_',REST_1_or_2,'_RL'],'Movement_RelativeRMS_mean.txt');
        movement_rel_RMS_RL = readmatrix(movement_rel_RMS_RL_path);
        is_excluded_due_movement(sub_idx,rest_idx) = ...
            any(movement_rel_RMS_LR > movement_thres) || any(movement_rel_RMS_RL > movement_thres);
    catch
        is_excluded_due_movement(sub_idx,rest_idx) = true;
    end
end

RL_processing_errors_subs = [103010;113417;116423;120010;121719;127226;130114;143830;169040;185038;189652;202820;204218;329844;385046;401422;462139;469961;644246;688569;723141;908860;943862;969476;971160];
is_RL_processing_errors = false(length(sub_ids),1);
for n_rl = 1:length(RL_processing_errors_subs)
    is_RL_processing_errors(sub_ids==RL_processing_errors_subs(n_rl)) = true;
end

%%
save results/HCP_timeseries_subject_exclude_info sub_ids is_excluded_due_movement does_have_MMSE is_cognitive_impaired MMSE_thres movement_thres is_RL_processing_errors
