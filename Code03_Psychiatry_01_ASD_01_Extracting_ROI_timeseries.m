% clear; clc;
current_path = pwd;

%%
data_ASD_path = 'G:\ABIDE1_func_preproc\Outputs\cpac\filt_noglobal\func_preproc';

siteList={'CALTECH','CMU','KKI','MAX_MUN','NYU','OLIN','OHSU','SDSU','SBL','STANFORD','TRINITY','UCLA','LEUVEN','UM','PITT','USM','YALE'};
TRlist=[2000 2000 2500 3000 2000 1500 2500 2000 2200 2000 2000 3000 1656 2000 1500 2000 2000];

label_cortical_file = 'atlas/MMP_cortex_resample.nii';
label_subcortical_file = 'atlas/CortexSubcortex_nifti_resample.nii';

label_cortex = niftiread(label_cortical_file);
label_subcortical = niftiread(label_subcortical_file);

label_subcortical_file_original = 'atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_nifti.nii';
label_subcortical_original = niftiread(label_subcortical_file_original);
%%
cd(data_ASD_path);
image_file_list = dir('*_func_preproc.nii.gz');
image_name_list = {image_file_list.name};
cd(current_path);

%%
roi_cortical_labels = unique(label_cortex(:));
roi_cortical_labels(roi_cortical_labels == 0) = []; % Remove background label if 0
roi_subcortical_labels = unique(label_subcortical(:));
roi_subcortical_labels(roi_subcortical_labels == 0) = []; % Remove background label if 0
roi_subcortical_labels_original = unique(label_subcortical_original(:));
roi_subcortical_labels_original(roi_subcortical_labels_original == 0) = []; % Remove background label if 0

time_series_denoised_filtered = cell(length(image_file_list),1);

for n_sub = 1:length(image_file_list)
    try
        tic
        fmri_data = niftiread([image_file_list(n_sub).folder, filesep, image_file_list(n_sub).name]);

        disp(n_sub);
        
        [dimX_fmri, dimY_fmri, dimZ_fmri, dimT_fmri] = size(fmri_data);
        fmri_reshaped = reshape(fmri_data, [], dimT_fmri);

        roi_time_series = zeros(718,dimT_fmri);
        for ii = 1:length(roi_cortical_labels)
            label = roi_cortical_labels(ii);
            roi_mask = label_cortex == label;
            voxel_indices = find(roi_mask);

            roi_time_series(ii,:) = mean(fmri_reshaped(voxel_indices, :));
        end

        for jj = 1:length(roi_subcortical_labels_original)
            if any(roi_subcortical_labels_original(jj) == roi_subcortical_labels)
                label = roi_subcortical_labels_original(jj);
                roi_mask = label_subcortical == label;
                voxel_indices = find(roi_mask);

                roi_time_series(jj+360,:) = mean(fmri_reshaped(voxel_indices, :));
            end
        end

        time_series_denoised_filtered{n_sub} = roi_time_series;
        toc
    catch
    end
end
%%
save results/ABIDE_timeseires_cortical_subcortical_CPAC_denoised_filtered time_series_denoised_filtered image_file_list siteList TRlist