clear; clc;
current_path = pwd;
load('results/HCP_timeseries_subject_exclude_info.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered_meta.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered.mat');

%%
is_sub_exclude = true;
if is_sub_exclude
    for nsub = 1:length(sub_ids)
        if does_have_MMSE(nsub) || is_cognitive_impaired(nsub) || is_RL_processing_errors(nsub)
            time_series_denoised_filtered(nsub,:) = {[],[],[],[]}; %#ok<SAGROW>
        else
            if is_excluded_due_movement(nsub,1)
                time_series_denoised_filtered(nsub,1:2) = {[],[]}; %#ok<SAGROW>
            end
            if is_excluded_due_movement(nsub,2)
                time_series_denoised_filtered(nsub,3:4) = {[],[]}; %#ok<SAGROW>
            end
        end
    end
end

remaining_sub_idx = false(length(sub_ids),1);
for nsub = 1:length(sub_ids)
    if ~isempty(time_series_denoised_filtered{nsub,1}) || ~isempty(time_series_denoised_filtered{nsub,2}) || ~isempty(time_series_denoised_filtered{nsub,3}) || ~isempty(time_series_denoised_filtered{nsub,4})
        remaining_sub_idx(nsub) = true;
    end
end

load secure_data/path_info;
gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
freesurfer_data_table = readtable(freesurfer_data_path);
for nrow = size(gene_data_table,1):-1:1
    if ~any(sub_ids==gene_data_table(nrow,'Subject').Variables)
        gene_data_table(nrow,:) = [];
    end
end
for nrow = size(behav_data_table,1):-1:1
    if ~any(sub_ids==behav_data_table(nrow,'Subject').Variables)
        behav_data_table(nrow,:) = [];
    end
end
gene_data_table = sortrows(gene_data_table, 'Subject');
behav_data_table = sortrows(behav_data_table, 'Subject');

ages = gene_data_table.Age_in_Yrs;
genders = behav_data_table.Gender;

num_female = sum(strcmp(genders(remaining_sub_idx),'F'));
mean_age = mean(ages(remaining_sub_idx));
std_age = std(ages(remaining_sub_idx));

fprintf('Total number of subjects: %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
    sum(remaining_sub_idx),num_female,mean_age,std_age);


%% Resting-state FC gradient
N=360;

corr_mats = zeros(4*size(time_series_denoised_filtered,1),N,N);
i_num = 0;
for ii = 1:(4*size(time_series_denoised_filtered,1))
    disp(ii);
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if ~isempty(time_series_denoised_filtered{nsub,nses})
        i_num = i_num + 1;
        corr_mat = corrcoef(time_series_denoised_filtered{nsub,nses}');
        corr_mats(i_num,:,:) = corr_mat(1:N,1:N);
    end
end

corr_mats = atanh(corr_mats(1:i_num,:,:));

try 
    label_idx_list;
catch
    label_idx_list=1:N;
end

group_corr_mat = tanh(squeeze(nanmean(corr_mats)));
group_corr_mat_original = tanh(squeeze(nanmean(corr_mats)));

[~,sort_idx] = sort(abs(group_corr_mat(:)),'descend');

group_corr_mat(sort_idx(length(sort_idx)*0.2:end)) = 0;


[coeff, score, latent, tsquared, explained, mu] = pca(group_corr_mat,'NumComponents',5);

RSFC_grad = score;
save results/RSFC_standard_grad RSFC_grad group_corr_mat group_corr_mat_original

%% Resting-state FC gradient (including subcortex)
N=718;

corr_mats = zeros(4*size(time_series_denoised_filtered,1),N,N);
i_num = 0;
for ii = 1:(4*size(time_series_denoised_filtered,1))
    disp(ii);
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if ~isempty(time_series_denoised_filtered{nsub,nses})
        i_num = i_num + 1;
        corr_mat = corrcoef(time_series_denoised_filtered{nsub,nses}');
        corr_mats(i_num,:,:) = corr_mat(1:N,1:N);
    end
end

corr_mats = atanh(corr_mats(1:i_num,:,:));

try 
    label_idx_list;
catch
    label_idx_list=1:N;
end

group_corr_mat = tanh(squeeze(nanmean(corr_mats)));
group_corr_mat_original = tanh(squeeze(nanmean(corr_mats)));

[~,sort_idx] = sort(abs(group_corr_mat(:)),'descend');

group_corr_mat(sort_idx(length(sort_idx)*0.2:end)) = 0;


[coeff, score, latent, tsquared, explained, mu] = pca(group_corr_mat,'NumComponents',5);

RSFC_grad = score;
save results/RSFC_standard_grad_all RSFC_grad group_corr_mat group_corr_mat_original

%% Global signal removed
N = 360;
corr_mats_global_removed = zeros(4*size(time_series_denoised_filtered,1),N,N);
i_num = 0;
for ii = 1:(4*size(time_series_denoised_filtered,1))
    disp(ii);
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if ~isempty(time_series_denoised_filtered{nsub,nses})
        i_num = i_num + 1;
        
        time_series = time_series_denoised_filtered{nsub,nses}';
        Gs = mean(time_series(:,1:N),2);
        
        X = [ones(length(Gs), 1), Gs];
        
        beta = X \ time_series;

        time_series_fitted = X * beta;

        time_series_residuals = time_series - time_series_fitted;
        
        corr_mat = corrcoef(time_series_residuals);
        corr_mats_global_removed(i_num,:,:) = corr_mat(1:N,1:N);
    end
end

corr_mats_global_removed = atanh(corr_mats_global_removed(1:i_num,:,:));
group_corr_mat_global_removed = tanh(squeeze(nanmean(corr_mats_global_removed)));

group_corr_mat_original = group_corr_mat_original(1:N,1:N);

save results/RSFC group_corr_mat_original group_corr_mat_global_removed

%% Resting-state FC gradient (subcortex)
load results/RSFC_standard_grad_all
load('DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_subExclude.mat', 'roi_exclude')

group_corr_mat_original(roi_exclude,:) = [];
group_corr_mat_original(:,roi_exclude) = [];

labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
labels_subcortical_data = labels_subcortical.cdata;
label_idx_list = unique(labels_subcortical_data);
label_idx_list(label_idx_list==0) = [];
label_exclude = label_idx_list(roi_exclude);
label_idx_list(roi_exclude) = [];
models = labels_subcortical.diminfo{1}.models;

subcortical_areas = {'hippocampus','thalamus','striatum','cerebellum','amygdala','brainstem'};
structure_model_idx_list = {[14,15],[20,21],[3,4,8,9,18,19],[10,11],[5,6],7};
for n_area = 1:length(subcortical_areas)
    structure_model_idx = structure_model_idx_list{n_area};
    label_select = [];
    for n_model = structure_model_idx
        label_select = [label_select;
            labels_subcortical.cdata(models{1,n_model}.start:(models{1,n_model}.start+models{1,n_model}.count-1))];
    end
    label_select = unique(label_select);
    target_roi_idx = zeros(1,length(label_select));
    target_exclude = false(1,length(label_select));
    for n_roi = 1:length(label_select)
        try
            target_roi_idx(n_roi) = find(label_idx_list == label_select(n_roi));
        catch
            target_exclude(n_roi) = true;
        end
    end
    label_select(target_exclude) = [];
    target_roi_idx(target_exclude) = [];

    C = group_corr_mat_original(target_roi_idx,:);
    C(:,target_roi_idx) = [];
    
%     gm = GradientMaps('kernel','na','approach','pca','alignment','','random_state',10,'verbose',true);
%     
%     gm = gm.fit(C);
%
%     coeff = gm.gradients{1};
    
    S = corrcoef(C');

    [coeff, score, latent, tsquared, explained, mu] = pca(S,'NumComponents',5);
    fprintf('%s -- Variance explained by PC1: %.3f, PC2:  %.3f \n',subcortical_areas{n_area},explained(1),explained(2))
    
    if n_area ~= 2
        score = -score;
    end
    

    gradient_cifti_data = zeros(size(labels_subcortical_data));
    for n_roi = 1:length(label_select)
        gradient_cifti_data(labels_subcortical_data==label_select(n_roi)) = score(n_roi,1);
    end
    fprintf('Plot info > min: %.5f, max: %.5f \n',min(score(:,1)), max(score(:,1)));
    plot_subcortex(gradient_cifti_data,['results/RSFC_',subcortical_areas{n_area},'_grad.jpg'],[],min(coeff(:,1)),max(coeff(:,1)),subcortical_areas{n_area});
    
    eval(['RSFC_grad_',subcortical_areas{n_area},'=score;']);
    eval(['target_roi_idx_',subcortical_areas{n_area},'=target_roi_idx;']);

end

save results/RSFC_subcortical_grad RSFC_grad_* subcortical_areas structure_model_idx_list target_roi_idx_*

%% Yeo 7 network partition
load results/RSFC
% Parameters
numClusters = 7;                     % Number of networks to identify
numROIs = size(group_corr_mat_global_removed, 1); % Number of Regions of Interest (ROIs)

% Output Files
networkAssignmentFile = 'ROI_Network_Assignments.mat';

% Network Labels based on Yeo et al. (2011)
networkLabels = {'Default Mode Network (DMN)', ...
                'Dorsal Attention Network (DAN)', ...
                'Ventral Attention Network (VAN)', ...
                'Limbic Network', ...
                'Frontoparietal Control Network', ...
                'Somatomotor Network', ...
                'Visual Network'};
            
% Check if the number of network labels matches the number of clusters
if length(networkLabels) ~= numClusters
    error('Number of network labels (%d) does not match the number of clusters (%d).', ...
        length(networkLabels), numClusters);
end

% Optional: Load ROI Spatial Coordinates for Visualization
% roiCoordinates should be a [numROIs x 3] matrix containing X, Y, Z coordinates
% Uncomment and modify the following line based on your data
% load('roi_coordinates.mat', 'roiCoordinates'); 

% Set random seed for reproducibility
rng(123);

%%% ------------------ Functional Connectivity Profile Preparation ------------------ %%
fprintf('Preparing functional connectivity profiles for clustering...\n');

% Each ROI's connectivity profile is represented by its row in the connectivity matrix
connectivityProfiles = group_corr_mat_global_removed;

% Optional: Normalize connectivity profiles (Z-score normalization across ROIs)
% This step ensures that each ROI's connectivity profile has zero mean and unit variance
connectivityProfiles_norm = zscore(connectivityProfiles, 0, 2);

% Handle any NaN or Inf values resulting from normalization
connectivityProfiles_norm(isnan(connectivityProfiles_norm)) = 0;
connectivityProfiles_norm(isinf(connectivityProfiles_norm)) = 0;

fprintf('Functional connectivity profiles prepared.\n');

%%% ------------------ Clustering for 7-Network Parcellation ------------------ %%
fprintf('Starting k-means clustering for 7-network parcellation...\n');

% Perform k-means clustering
% Number of clusters is set to 7 as per Yeo et al. (2011)
opts = statset('Display','final','MaxIter',1000);
[clusterLabels, clusterCenters] = kmeans(connectivityProfiles_norm, numClusters, ...
    'Replicates', 10, 'Options', opts);

fprintf('Clustering completed.\n');
%%% ------------------ Save visulization ------------------ %%
N = 360;
label_idx_list = 1:N;
%%% save files
for k = 1:numClusters
    fprintf('Network %d created.\n', k);
    
    labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

    grad1_data = zeros(size(labels.cdata));
    for n_reg = 1:N
        grad1_data(labels.cdata==label_idx_list(n_reg)) = clusterCenters(k,n_reg);
    end
    lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
    cifti_write(lag1, ['results/Yeo_7Network_',num2str(k),'.dscalar.nii']);
end

%%% ------------------ Network Assignment and Labeling ------------------ %%
Yeo_7network_DMN1_vector = clusterCenters(3,:)';
Yeo_7network_SM_vector = clusterCenters(1,:)';
Yeo_7network_CEN_vector = clusterCenters(5,:)';
Yeo_7network_Vis_vector = clusterCenters(2,:)';
Yeo_7network_DAN_vector = clusterCenters(4,:)';
Yeo_7network_VAN_vector = clusterCenters(7,:)';
Yeo_7network_DMN2_vector = clusterCenters(6,:)';

Yeo_7network_all_vector = clusterCenters';

save results/Yeo_7networks Yeo_7network_*_vector

%%% ------------------ Clustering for 17-Network Parcellation ------------------ %%
fprintf('Starting k-means clustering for 17-network parcellation...\n');
numClusters = 17;
% Perform k-means clustering
% Number of clusters is set to 7 as per Yeo et al. (2011)
opts = statset('Display','final','MaxIter',1000);
[~, clusterCenters] = kmeans(connectivityProfiles_norm, numClusters, ...
    'Replicates', 10, 'Options', opts);

fprintf('Clustering completed.\n');
%%% ------------------ Save visulization ------------------ %%
N = 360;
label_idx_list = 1:N;
%%% save files
for k = 1:numClusters
    fprintf('Network %d created.\n', k);
    
    labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

    grad1_data = zeros(size(labels.cdata));
    for n_reg = 1:N
        grad1_data(labels.cdata==label_idx_list(n_reg)) = clusterCenters(k,n_reg);
    end
    lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
    cifti_write(lag1, ['results/Yeo_17Network_',num2str(k),'.dscalar.nii']);
end

%%% ------------------ Network Assignment and Labeling ------------------ %%

Yeo_17network_all_vector = clusterCenters';

save results/Yeo_17networks Yeo_17network_*_vector
%% PCC (31pv) seed FC (Default mode)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

precu_seed = mean(group_corr_mat_global_removed(1:N,[35,35+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

grad1_data = zeros(size(labels.cdata));
for n_reg = 1:N
    grad1_data(labels.cdata==label_idx_list(n_reg)) = precu_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
cifti_write(lag1, ['results/RSFC_Task_negative.dscalar.nii']);

Task_negative_vector = precu_seed;
save results/Task_negative_vector Task_negative_vector
    
%% SMG seed FC (Fronto-parietal)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_global_removed(1:N,[148,148+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
FPN_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_PFN.dscalar.nii']);

save results/FPN_vector FPN_vector

%% Somatosensory seed FC
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_global_removed(1:N,[52,52+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
SMLV_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_SMLV.dscalar.nii']);

save results/SMLV_vector SMLV_vector

%% aIns seed FC (Salience or Cingulo-Opercular)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_global_removed(1:N,[111]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
Salience_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_Salience.dscalar.nii']);

save results/Right_aIns_seed_vector Salience_vector

%% FEF seed FC (Dorsal attention)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_global_removed(1:N,[10,10+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
DAN_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_DAN.dscalar.nii']);

save results/DAN_vector DAN_vector

%% TPJ seed FC (ventral attention)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_global_removed(1:N,[139,139+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
VAN_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_VAN.dscalar.nii']);

save results/VAN_vector VAN_vector

%% Language seed FC
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_global_removed(1:N,[205,255,311]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
LAN_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_LAN.dscalar.nii']);

save results/LAN_vector LAN_vector

%% Lateralized index
N=360;

frame_dt = 0.72;
time_window_length = 30;
time_window_interval = 10;

MLI_all = zeros(4*size(time_series_denoised_filtered,1),N);
i_num = 0;
for ii = 1:(4*size(time_series_denoised_filtered,1))
    disp(ii);
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if ~isempty(time_series_denoised_filtered{nsub,nses})
        i_num = i_num + 1;
        
        time_series_matrix = time_series_denoised_filtered{nsub,nses};
        time_series_matrix = time_series_matrix(1:N,:);
        
        GS_R = mean(time_series_matrix(1:N/2,:));
        GS_L = mean(time_series_matrix(N/2+1:end,:));
        
        DLI = zeros(N,1);
    
        current_frame = 1;
        current_time = current_frame*frame_dt;
        accum_num = 0;
        while current_time+time_window_length < frame_dt*size(time_series_matrix,2)
            accum_num = accum_num + 1;

            current_GS_L = GS_L(current_frame:round(current_frame+time_window_length/frame_dt));
            current_GS_R = GS_R(current_frame:round(current_frame+time_window_length/frame_dt));

            current_time_series_matrix = time_series_matrix(:,current_frame:round(current_frame+time_window_length/frame_dt));

            r_L = corrcoef([current_GS_L',current_time_series_matrix']);
            r_R = corrcoef([current_GS_R',current_time_series_matrix']);

            DLI(1:N,accum_num) = (r_L(1,2:end) - r_R(1,2:end))';

            current_frame = current_frame + round(time_window_interval/frame_dt);
            current_time = current_frame*frame_dt;
        end

        MLI = mean(DLI,2);
        
        MLI_all(i_num,:) = MLI;
    end
end

MLI_all = MLI_all(1:i_num,:)';
MLI = mean(MLI_all,2);

save results/Lateralized_index MLI_all MLI

labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

grad1_data = zeros(size(labels.cdata));
for n_reg = 1:N
    grad1_data(labels.cdata==label_idx_list(n_reg)) = MLI(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
cifti_write(lag1, ['results/MLI.dscalar.nii']);

%% Time lag structure
% Define constants
num_nodes = 360;  % Define ROI count
max_sessions = 4; % Maximum number of sessions per subject
% Determine the number of subjects and sessions from the data
[num_subs, ~] = size(time_series_denoised_filtered);

% Set analysis parameters
lags = -6:6;        % Range of TR shifts
tr = 0.72;          % Sampling interval in seconds
lag_lim = 4;        % Lag limit (in seconds)


% Preallocate group matrices to accommodate all possible sessions
% Assuming maximum total sessions is n_sub * max_sessions
% Initialize with NaNs to handle missing data gracefully
grp_lags = single(nan(num_nodes, num_nodes, num_subs * max_sessions));
grp_AWTD = single(nan(num_nodes, num_nodes, num_subs * max_sessions));
grp_ZL = single(nan(num_nodes, num_nodes, num_subs * max_sessions));
grp_peak = single(nan(num_nodes, num_nodes, num_subs * max_sessions));

% Initialize NaN counters for group matrices
grp_lags_nans = single(zeros(num_nodes, num_nodes, num_subs * max_sessions));
grp_AWTD_nans = single(zeros(num_nodes, num_nodes, num_subs * max_sessions));
grp_ZL_nans = single(zeros(num_nodes, num_nodes, num_subs * max_sessions));
grp_peak_nans = single(zeros(num_nodes, num_nodes, num_subs * max_sessions));

% Initialize session counter
session_counter = 0;

% Loop over subjects
for n_sub = 1:num_subs
    % Loop over sessions for the current subject
    for n_ses = 1:max_sessions
        tic
        disp(['Processing Subject ' num2str(n_sub) ', Session ' num2str(n_ses)]);
        
        % Load the current session's BOLD data
        BOLD = time_series_denoised_filtered{n_sub, n_ses};
        
        % Check if the BOLD data exists and is non-empty
        if isempty(BOLD) || size(BOLD, 1) == 0 || size(BOLD, 2) == 0
            disp(['  Skipping Session ' num2str(n_ses) ' for Subject ' num2str(n_sub) ' due to empty data']);
            continue;  % Skip to the next session
        end
        
        BOLD = BOLD(1:num_nodes,:);
        
        disp(['  BOLD size: ', mat2str(size(BOLD))]);  % BOLD size check
        
        % De-mean the time series (for covariance calculation)
        run_mean = nanmean(BOLD, 2);
        BOLD = bsxfun(@minus, BOLD, run_mean);
        
        % Compute the lagged covariance (TD matrices)
        Cov = lagged_cov(BOLD', BOLD', max(lags));  % Cross-covariance calculation across ROIs
        disp(['  Covariance matrix size: ', mat2str(size(Cov))]);  % Covariance matrix size check
        
        % Normalize pairwise cross-covariance functions based on the entire run
        for k = 1:numel(lags)
            Cov(:,:,k) = Cov(:,:,k) / (size(BOLD, 1) - abs(lags(k)));  % Normalize by time axis
        end
        
        % Parabolic interpolation to get peak lag/correlation
        [pl, pc] = parabolic_interp(Cov, tr);
        pl(abs(pl) > lag_lim) = nan;  % Exclude long lags
        disp(['  pl size: ', mat2str(size(pl))]);  % Check size
        disp(['  pc size: ', mat2str(size(pc))]);  % Check size
        
        % Get zero-lag correlation
        temp = Cov(:,:,lags == 0);  % Zero-lag correlation
        d = zeros(size(temp));
        d(logical(eye(length(temp)))) = sqrt(diag(temp));  % Diagonal normalization
        temp = d^(-1) * temp * d^(-1);
        temp = atanh(temp);  % Fisher z-transform
        temp(isnan(pl)) = nan;  % Exclude NaNs from pl
        
        % Compute AWTD (Amplitude-Weighted Time Delay) matrix
        AWTD = pl .* pc;  % Weighting the time delays (pl) by peak correlation (pc)
        
        % Increment session counter
        session_counter = session_counter + 1;
        
        % Store computed matrices into the group matrices
        grp_lags(:,:,session_counter) = single(pl);
        grp_peak(:,:,session_counter) = single(pc);
        grp_ZL(:,:,session_counter) = single(temp);
        grp_AWTD(:,:,session_counter) = single(AWTD);
        
        % Accumulate NaNs for this session
        grp_lags_nans(:,:,session_counter) = single(isnan(pl));
        grp_peak_nans(:,:,session_counter) = single(isnan(pc));
        grp_ZL_nans(:,:,session_counter) = single(isnan(temp));
        grp_AWTD_nans(:,:,session_counter) = single(isnan(AWTD));
        
        toc
    end
end

% Trim group matrices to include only processed sessions
grp_lags = grp_lags(:,:,1:session_counter);
grp_peak = grp_peak(:,:,1:session_counter);
grp_ZL = grp_ZL(:,:,1:session_counter);
grp_AWTD = grp_AWTD(:,:,1:session_counter);

grp_lags_nans = grp_lags_nans(:,:,1:session_counter);
grp_peak_nans = grp_peak_nans(:,:,1:session_counter);
grp_ZL_nans = grp_ZL_nans(:,:,1:session_counter);
grp_AWTD_nans = grp_AWTD_nans(:,:,1:session_counter);

% Compute group averages across all valid sessions
grp_lags_mean = nanmean(grp_lags, 3);  % Mean across sessions for lags
grp_peak_mean = nanmean(grp_peak, 3);  % Mean across sessions for peak correlation
grp_ZL_mean = tanh(nanmean(grp_ZL, 3));  % Mean across sessions for zero-lag correlation, followed by un-Fisher z-transform
grp_AWTD_mean = nanmean(grp_AWTD, 3);  % Mean across sessions for AWTD

% Visualize results
figure; imagesc(grp_lags_mean, [-0.75, 0.75]); colorbar;
title('Group-Level Time Lag Matrix');

figure; imagesc(grp_peak_mean, [-0.7, 0.7]); colorbar;
title('Group-Level Peak Correlation Matrix');

figure; imagesc(grp_ZL_mean, [-0.7, 0.7]); colorbar;
title('Group-Level Zero-Lag Correlation Matrix');

figure; imagesc(grp_AWTD_mean, [-0.7, 0.7]); colorbar;
title('Group-Level Amplitude-Weighted Time Delay Matrix (AWTD)');

save results/time_lags grp_lags_mean grp_peak_mean grp_ZL_mean grp_AWTD_mean

%%
% Calculate Group-Level Coactivation Patterns (CAPs) from fMRI Data
% Author: Youngjo Song
% Date: 2024-09-30

%%% Parameters (Adjust as Needed)
numClusters = 8;                  % Number of clusters for k-means
thresholdPercentile = 85;         % Percentile for frame selection
seedIndices = [35, 35+180];     % Indices of the PCC (31pv) seed regions
numSubjects = 1096;               % Total number of subjects
numSessions = 4;                  % Number of sessions per subject
maxFramesPerSubject = 200;        % Maximum frames to select per subject-session to manage memory
CAPsSaveFile = 'results/Group_CAPs.mat';  % File to save the resulting CAPs
selectedFrameAccumulator = [];    % Initialize accumulator for selected frames
missingDataList = {};             % Initialize cell array to log missing data
N = 360;

% Optionally, set random seed for reproducibility
rng(1);

%%% Step 1: Data Preprocessing and Frame Selection

fprintf('Starting data preprocessing and frame selection...\n');

for n_sub = 1:numSubjects
    for n_session = 1:numSessions
        % Display progress every 100 subjects at the start of their first session
        if mod(n_sub, 100) == 0 && n_session == 1
            fprintf('Processing Subject %d/%d...\n', n_sub, numSubjects);
        end
        
        % Extract fMRI data for current subject and session
        % Check if the cell is empty to handle missing data
        currentData = time_series_denoised_filtered{n_sub, n_session};
        if isempty(currentData)
            warning('Subject %d, Session %d: Data is missing. Skipping...', n_sub, n_session);
            missingDataList{end+1, 1} = n_sub;       % Log subject number
            missingDataList{end, 2} = n_session;     % Log session number
            continue; % Skip to next session
        end
        
        % Validate that currentData is a numeric matrix
        if ~isnumeric(currentData) || ~ismatrix(currentData)
            warning('Subject %d, Session %d: Data is not a valid numeric matrix. Skipping...', n_sub, n_session);
            missingDataList{end+1, 1} = n_sub;
            missingDataList{end, 2} = n_session;
            continue; % Skip to next session
        end
        
        currentData(361:end,:) = [];
        
        % Step 1.1: Data Preprocessing
        % Demean: Subtract the mean across time for each ROI
        fMRI_demeaned = currentData - mean(currentData, 2);
        
        % Normalize: Divide by the standard deviation across time for each ROI
        fMRI_std = std(fMRI_demeaned, 0, 2);
        fMRI_std(fMRI_std == 0) = 1; % Prevent division by zero
        fMRI_normalized = fMRI_demeaned ./ fMRI_std;
        
        % Handle any NaN or Inf values resulting from division
        fMRI_normalized(isnan(fMRI_normalized)) = 0;
        fMRI_normalized(isinf(fMRI_normalized)) = 0;
        
        % Step 1.2: Seed Selection
        % Extract the BOLD signal from the precuneus seed regions
        seedSignals = fMRI_normalized(seedIndices, :);
        
        % Average the seed signals if multiple seeds are provided
        seedAverage = mean(seedSignals, 1);
        
        % Step 1.3: Frame Selection
        % Determine the threshold based on the specified percentile of the seed signal
        threshold = prctile(seedAverage, thresholdPercentile);
        
        % Identify time points where the seed signal exceeds the threshold
        selectedTimePoints = find(seedAverage > threshold);
        
        % Optional: Limit the number of frames per subject-session to manage memory
        if length(selectedTimePoints) > maxFramesPerSubject
            selectedTimePoints = selectedTimePoints(randperm(length(selectedTimePoints), maxFramesPerSubject));
        end
        
        % Check if any frames are selected
        if isempty(selectedTimePoints)
            warning('Subject %d, Session %d: No frames exceed the threshold. Skipping...', n_sub, n_session);
            missingDataList{end+1, 1} = n_sub;
            missingDataList{end, 2} = n_session;
            continue; % Skip to next session if no frames are selected
        end
        
        % Extract the selected frames
        selectedFrames = fMRI_normalized(:, selectedTimePoints)'; % [Selected Frames x ROIs]
        
        % Accumulate selected frames
        selectedFrameAccumulator = [selectedFrameAccumulator; selectedFrames];
        
        % Clear variables to free memory
        clear currentData fMRI_demeaned fMRI_normalized seedSignals seedAverage threshold selectedTimePoints selectedFrames;
    end
end

fprintf('Data preprocessing and frame selection completed.\n');
fprintf('Total selected frames: %d\n', size(selectedFrameAccumulator, 1));

%%% Step 2: Clustering for Group-Level CAPs

fprintf('Starting k-means clustering on accumulated frames...\n');

% To manage memory and computational load, consider using a subset or
% dimensionality reduction (e.g., PCA) before clustering.

% Optional Step 2.1: Dimensionality Reduction with PCA
% Perform PCA without limiting the number of components to get all variances
[coeff_full, score_full, latent_full] = pca(selectedFrameAccumulator, 'Centered', true, 'NumComponents', size(selectedFrameAccumulator,2));

% Calculate cumulative variance
cumulativeVariance = cumsum(latent_full) / sum(latent_full);

% Find the number of components needed to retain at least 99.9% variance
numPCAComponents = find(cumulativeVariance >= 0.999, 1, 'first');

fprintf('Number of PCA components to retain 99.9% variance: %d\n', numPCAComponents);

[coeff, score, ~] = pca(selectedFrameAccumulator, 'Centered', true, 'NumComponents', numPCAComponents);


% Step 2.2: Perform k-means clustering on PCA-reduced data
opts = statset('Display','final','MaxIter',1000);
[clusterLabels, clusterCenters] = kmeans(score, numClusters, 'Replicates', 5, 'Options', opts);

fprintf('K-means clustering completed.\n');

%%% Step 3: CAP Creation

fprintf('Creating group-level CAPs...\n');

% Initialize matrix to hold CAPs
CAPs = zeros(numClusters, 360);

% Reconstruct cluster centers to original ROI space using PCA coefficients
reconstructedCenters = clusterCenters * coeff';

label_idx_list = 1:N;

for k = 1:numClusters
    % Each cluster center is already the mean in PCA space; reconstruct to ROI space
    CAPs(k, :) = reconstructedCenters(k, :);
    fprintf('CAP %d created.\n', k);
    
    labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

    grad1_data = zeros(size(labels.cdata));
    for n_reg = 1:N
        grad1_data(labels.cdata==label_idx_list(n_reg)) = CAPs(k,n_reg);
    end
    lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
    cifti_write(lag1, ['results/CAP_PCC_seed_',num2str(k),'.dscalar.nii']);
end

%%% Step 4: Saving and Visualization

% Save the group-level CAPs to a .mat file
save(CAPsSaveFile, 'CAPs');

fprintf('Group-level CAPs saved to %s.\n', CAPsSaveFile);

%%% End of Script
fprintf('Group-level CAP calculation completed successfully.\n');

%% QPP
% MATLAB Script to Estimate Quasi-Periodic Patterns (QPPs) Using Template Autoregressive Matching
% Based on Majeed et al.18 Algorithm
% Author: Youngjo Song
% Date: 2024/10/12

% Parameters
% num_subjects = 100; %size(time_series_denoised_filtered, 1); 
target_num_subjects = 100;
num_sessions = size(time_series_denoised_filtered, 2); % 4
num_rois = 360;
time_points_per_session = 1200;

% Initialize a list to store all concatenated data
% Preallocate for efficiency (assuming all cells are non-empty; adjust if necessary)
% Count non-empty cells
num_valid_cells = 0;
num_valid_subjects = 0;
for subj = 1:size(time_series_denoised_filtered,1)
    sub_valid_check = false;
    for sess = 1:num_sessions
        if ~isempty(time_series_denoised_filtered{subj, sess})
            num_valid_cells = num_valid_cells + 1;
            sub_valid_check = true;
        end
    end
    if sub_valid_check
        num_valid_subjects = num_valid_subjects + 1;
    end
    if num_valid_subjects == target_num_subjects
        break
    end
end

num_subjects = subj;

% Preallocate concatenated data matrix
% Size: 360 x (1200 * num_valid_cells)
group_data = zeros(num_rois, time_points_per_session * num_valid_cells);

% Concatenate data
current_index = 1;
for subj = 1:num_subjects
    disp(['Concatanating sub-',num2str(subj)]);
    for sess = 1:num_sessions
        if ~isempty(time_series_denoised_filtered{subj, sess})
            fmri_data = time_series_denoised_filtered{subj, sess}(1:num_rois, :); % 360x1200
            % Compute the global signal (mean across ROIs for each time point)
            global_signal = mean(fmri_data, 1)';

            % Center the global signal
            global_signal = global_signal - mean(global_signal);

            % Construct the design matrix with a constant term and the global signal
            X = [ones(length(global_signal), 1), global_signal];

            % Initialize the output data
            fmri_data_gsr = zeros(size(fmri_data));

            % Regress out the global signal from each ROI's time series
            for roi = 1:size(fmri_data, 1)
                y = fmri_data(roi, :)';
                beta = X \ y; % Ordinary Least Squares regression
                y_pred = X * beta;
                fmri_data_gsr(roi, :) = normalize((y - y_pred)');
            end

            group_data(:, current_index:(current_index + time_points_per_session -1)) = fmri_data_gsr;
%             group_data(:, current_index:(current_index + time_points_per_session -1)) = fmri_data;
            current_index = current_index + time_points_per_session;
        end
    end
end
clear time_series_denoised_filtered

% Trim any unused preallocated space (if any)
group_data = group_data(:, 1:(current_index -1));

% Transpose for easier windowing (Time x ROI)
group_data = group_data'; % Size: T x 360, where T = total time points

% Total time points after concatenation
T = size(group_data, 1);

%%% Define QPP Estimation Parameters
window_length = 30;        % Number of TRs in each window
initial_corr_threshold = 0.1;      % Correlation threshold to identify similar segments
final_corr_threshold = 0.2;      % Correlation threshold to identify similar segments
num_runs = 10;             % Number of independent runs
max_iterations = 50;       % Maximum iterations per run
convergence_threshold = 0.9999; % Threshold to determine convergence

%%% Initialize Variables to Store Run Results
templates = cell(num_runs, 1);              % To store templates from each run
avg_correlations = zeros(num_runs, 1);      % To store average correlations per run
corr_time_series_runs = cell(num_runs, 1);  % To store correlation time series per run

%%% QPP Estimation via Template Autoregressive Matching
fprintf('Starting QPP estimation with %d runs...\n', num_runs);
for run = 1:num_runs
    fprintf('Run %d/%d\n', run, num_runs);
    
    % Initialize with a random window
    rand_start = randi([1, T-window_length+1]);
    template_window = group_data(rand_start:(rand_start + window_length -1), :); % 30x360
    template_vector = template_window(:); % 10800x1 vector
    
    % Normalize the template vector
    template_vector = (template_vector - mean(template_vector)) / std(template_vector);
    
    corr_coeffs = zeros(T - window_length +1, 1);
    
    % Initialize variables for convergence
    converged = false;
    iteration = 0;
    
    while ~converged
        tic
        iteration = iteration +1;
        if iteration < 3
            current_corr_threshold = initial_corr_threshold;
        else
            current_corr_threshold = final_corr_threshold;
        end
        
        % Initialize a vector to store correlation coefficients
        corr_coeffs(:) = 0;
        
        % Compute correlation between template and all sliding windows
        fprintf('  Iteration %d: Computing correlations...\n', iteration);
        parfor s = 1:(T - window_length +1)
            current_window = group_data(s:(s + window_length -1), :); % 30x360
            current_vector = current_window(:); % 10800x1
            
            % Normalize the current window vector
            current_vector = (current_vector - mean(current_vector)) / std(current_vector);
            
            % Compute Pearson correlation
            r = corr(template_vector, current_vector);
            corr_coeffs(s) = r;
            
            % Optional: Progress indicator every 100,000 windows
%             if mod(s, 100000) == 0
%                 fprintf('    Processed %d/%d windows...\n', s, T - window_length +1);
%             end
        end
        
        % Identify windows with correlation > threshold
        [~, similar_indices] = findpeaks(corr_coeffs, 'MinPeakHeight', current_corr_threshold);
%         similar_indices = find(corr_coeffs > corr_threshold);
        num_similar = length(similar_indices);
        fprintf('    Found %d windows with local maxima with r > %.2f\n', num_similar, current_corr_threshold);
        
        if num_similar ==0
            warning('    No windows found with correlation above threshold. Stopping iterations.');
            break;
        end
        
        % Extract all similar windows and average them to update the template
        similar_windows = zeros(num_similar, window_length * num_rois); % Preallocate
        for idx =1:num_similar
            s = similar_indices(idx);
            window = group_data(s:(s + window_length -1), :); % 30x360
            similar_windows(idx, :) = window(:)';
        end
        
        % Compute the new template by averaging similar windows
        new_template_vector = mean(similar_windows, 1)'; % 10800x1
                
        % Check for convergence (change in template below threshold)
%         delta = norm(new_template_vector - template_vector);
        cc = corrcoef(new_template_vector,template_vector);
        fprintf('    Iteration %d: Similarity in template = %.6f\n', iteration, cc(2));
        if cc(2) > convergence_threshold
            fprintf('    Convergence reached.\n');
            converged = true;
        end
        
        % Update the template
        template_vector = new_template_vector;
        toc
    end
    
    % Store the converged template
    final_template = reshape(template_vector, [window_length, num_rois])'; % 360x30
    templates{run} = final_template;
    
    % Compute the average correlation for this run
    avg_correlations(run) = mean(corr_coeffs);
%     corr_time_series_runs{run} = corr_coeffs;
    
    fprintf('  Run %d completed. Average correlation: %.4f\n', run, avg_correlations(run));
end

%%% Select the Best Template Based on Average Correlation
[~, best_run] = max(avg_correlations);
final_qpp_template = templates{best_run}; % 360x30
% final_corr_time_series = corr_time_series_runs{best_run}; % (T -29)x1

fprintf('QPP estimation completed. Best run: %d with average correlation: %.4f\n', best_run, avg_correlations(best_run));

%%% Visualization of the Final QPP Template and Correlation Time Series

% Compute the global signal of the final QPP template
% Global signal: mean across ROIs for each time point in the window
final_qpp_global = mean(final_qpp_template, 1); % 1x30

% Plot the final QPP global signal
figure;
plot(final_qpp_global, 'LineWidth', 2);
xlabel('Time Points (TRs)');
ylabel('Normalized Global Signal');
title('Final QPP Global Signal');
grid on;


disp('QPP estimation using the template autoregressive matching algorithm completed successfully.');

save results/QPP_100_GSR final_qpp_template avg_correlations templates

%% display QPP
N_cortex = 360;
cd(current_path);
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

atlasimage=niftiread('atlas/MMP_in_MNI_symmetrical_LR_diff.nii');
atlasinfo=niftiinfo('atlas/MMP_in_MNI_symmetrical_LR_diff.nii');

atlasimage_temp = atlasimage;
for n_roi = 1:(N_cortex/2)
    atlasimage(atlasimage_temp==n_roi) = n_roi + 180;
end
for n_roi = ((N_cortex/2)+1):N_cortex
    atlasimage(atlasimage_temp==n_roi) = n_roi - 180;
end

try
    label_idx_list;
catch
    label_idx_list = 1:N;
end

mkdir('results/QPP_GSR');
save_dir = [pwd filesep 'results/QPP_GSR'];

frame_dt = 0.72;
N = 360;

cd(save_dir);
for frame = 1:size(final_qpp_template,2)
    eigenstate = zeros(size(atlasimage));
    for roi=1:N 
        eigenstate(atlasimage==label_idx_list(roi)) = final_qpp_template(roi,frame);
    end

    atlasinfo.Datatype='double';
    niftiwrite(eigenstate,['QPP_template_', num2str(frame,'%04d')],atlasinfo);

    fh = conn_mesh_display(['QPP_template_', num2str(frame,'%04d'), '.nii']);
    fh('colormap','bluewhitered');
    fh('colorbar','on');
    fh('colorbar','rescale',[min(final_qpp_template,[],'all'),max(final_qpp_template,[],'all')]);
    fh('print',4,['QPP_template_', num2str(frame,'%04d'),'.jpg'],'-r150','-nogui') 

    close;
end


% Define the folder where your images are stored
imageFolder = pwd;
imageFiles = dir(fullfile(imageFolder, '*.jpg')); % Adjust the pattern if needed

niiFiles = dir(fullfile(imageFolder, '*.nii'));
for ii = 1:length(niiFiles)
    delete([niiFiles(ii).folder,filesep,niiFiles(ii).name]);
end

% Create a VideoWriter object
outputVideo = VideoWriter(fullfile(imageFolder, 'outputVideo.avi'));
outputVideo.FrameRate = 10; % Set frame rate

% Open the video writer
open(outputVideo);

% Loop through each image, read it, and write it to the video
for i = 1:length(imageFiles)
    img = imread(fullfile(imageFolder, imageFiles(i).name));
    if i == 1
        desiredSize = [size(img,1), size(img,2)];
        resizedImg = img;
    else
        resizedImg = imresize(img, desiredSize); % Resize the image
    end

    % Calculate current time
    currentTime = 0 + (i - 1) * frame_dt;
    timeText = sprintf('t=%.1f', currentTime); % Format text

    % Add text to the image
    position = [20, 20]; % Position of the text (x, y)
    fontSize = 50; % Font size
    annotatedImg = insertText(double(resizedImg)/255, position, timeText, 'FontSize', fontSize, 'TextColor', 'black', 'BoxOpacity', 0);

    % Write the annotated image to the video
    writeVideo(outputVideo, annotatedImg);
end

% Close the video writer
close(outputVideo);

cd(current_path);
