clear; clc;
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered_meta.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered.mat');
%% Resting-state FC gradient
N=360;
% N=718;

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

RSFC_grad = coeff;
save results/RSFC_standard_grad RSFC_grad group_corr_mat group_corr_mat_original

%% Global signal removed
corr_mats_global_removed = zeros(4*size(time_series_denoised_filtered,1),N,N);
i_num = 0;
for ii = 1:(4*size(time_series_denoised_filtered,1))
    disp(ii);
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if ~isempty(time_series_denoised_filtered{nsub,nses})
        i_num = i_num + 1;
        
        time_series = time_series_denoised_filtered{nsub,nses}';
        Gs = mean(time_series,2);
        
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
%% Precuneuous seed FC (Default mode)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

precu_seed = mean(group_corr_mat_global_removed(1:N,[161,161+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

grad1_data = zeros(size(labels.cdata));
for n_reg = 1:N
    grad1_data(labels.cdata==label_idx_list(n_reg)) = precu_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
cifti_write(lag1, ['results/RSFC_Task_positive.dscalar.nii']);

DMN_vector = precu_seed;
save results/DMN_vector DMN_vector
    
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

area_seed = mean(group_corr_mat_original(1:N,[52,52+180]),2);
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

area_seed = mean(group_corr_mat_original(1:N,[111,111+180]),2);
% precu_seed(abs(precu_seed)<0.2) = 0;

FPN_grad = zeros(size(labels.cdata));
Salience_vector = area_seed;
for n_reg = 1:N
    FPN_grad(labels.cdata==label_idx_list(n_reg)) = area_seed(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, FPN_grad, 'dscalar');
cifti_write(lag1, ['results/RSFC_Salience.dscalar.nii']);

save results/Salience_vector Salience_vector

%% FEF seed FC (Dorsal attention)
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

area_seed = mean(group_corr_mat_original(1:N,[10,10+180]),2);
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

area_seed = mean(group_corr_mat_original(1:N,[139,139+180]),2);
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

area_seed = mean(group_corr_mat_original(1:N,[205,255,311]),2);
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
        
        GS_L = mean(time_series_matrix(1:N/2,:));
        GS_R = mean(time_series_matrix(N/2+1:end,:));
        
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
    grad1_data(labels.cdata==label_idx_list(n_reg)) = MSI(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
cifti_write(lag1, ['results/MLI.dscalar.nii']);