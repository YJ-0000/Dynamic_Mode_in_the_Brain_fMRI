clear;clc;
current_path = pwd;
%%
labels_ref = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
labels_cortical = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
label_cortical_data = round(labels_cortical.cdata);
labels_subcortical = cifti_read('atlas/Cortical_Subcortical.dscalar.nii');

try
    label_subcortical_data = round(labels_subcortical.cdata);
    label_data = [label_cortical_data; label_subcortical_data(length(label_cortical_data)+1:end)];
catch
    label_data = label_cortical_data;
end

HCP_denoised_path = 'G:\HCP_3T_denoised';
HCP_preprocessed_path = 'G:\HCP_3T_preproc';
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
save HCP_timeseries_cortical_subcortical_extracted_meta sub_ids folder_denoised_list N len_time label_idx_list
save('HCP_timeseries_cortical_subcortical_extracted.mat', 'time_series_denoised', '-v7.3')