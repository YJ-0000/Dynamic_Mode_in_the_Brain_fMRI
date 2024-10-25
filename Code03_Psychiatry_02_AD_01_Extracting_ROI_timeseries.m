clear; clc;
current_path = pwd;

data_path = 'D:\ADNI';

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

%% structural normalization
cd(data_path);
data_dir = dir;
data_dir(~[data_dir.isdir]) = [];
data_dir(1:2) = [];

for ndir = 1:length(data_dir)
    cd([data_dir(ndir).folder,filesep,data_dir(ndir).name]);
    sub_dir = dir('*_S_*');
    for nsub = 1:length(sub_dir)
        fprintf(['sub >> ',sub_dir(nsub).name,'.\n']);
        cd([sub_dir(nsub).folder,filesep,sub_dir(nsub).name]);
        cd('anat');
        deform_file = dir('y_*.nii');
        wm_file = dir('c2*.nii');
        csf_file = dir('c3*.nii');
        
        matlabbatch = cell(1);
        matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[deform_file(1).folder,filesep,deform_file(1).name]};
        matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[wm_file(1).folder,filesep,wm_file(1).name]
                                                                    [csf_file(1).folder,filesep,csf_file(1).name]};
        matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                  78 76 85];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
        
        spm_jobman('run',matlabbatch);
    end
end

%% ROI extraction
cd(current_path);
label_cortical_file = 'atlas/MMP_cortex_resample_AD.nii';
label_subcortical_file = 'atlas/CortexSubcortex_nifti_resample_AD.nii';

label_cortex = niftiread(label_cortical_file);
label_subcortical = niftiread(label_subcortical_file);

label_subcortical_file_original = 'atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR_nifti.nii';
label_subcortical_original = niftiread(label_subcortical_file_original);

roi_cortical_labels = unique(label_cortex(:));
roi_cortical_labels(roi_cortical_labels == 0) = []; % Remove background label if 0
roi_subcortical_labels = unique(label_subcortical(:));
roi_subcortical_labels(roi_subcortical_labels == 0) = []; % Remove background label if 0
roi_subcortical_labels_original = unique(label_subcortical_original(:));
roi_subcortical_labels_original(roi_subcortical_labels_original == 0) = []; % Remove background label if 0

AD_timeseries_extracted = cell(1,length(data_dir));
sub_list = cell(1,length(data_dir));
for ndir = 1:length(data_dir)
    cd([data_dir(ndir).folder,filesep,data_dir(ndir).name]);
    sub_dir = dir('*_S_*');
    temp_timeseries_extracted = cell(length(sub_dir),1); % ROI, movement, WM, CSF
    temp_sub_list = cell(length(sub_dir),1);
    for nsub = 1:length(sub_dir)
        fprintf(['sub >> ',sub_dir(nsub).name,'.\n']);
        cd([sub_dir(nsub).folder,filesep,sub_dir(nsub).name]);
        cd('anat');
        wm_file = dir('wc2*.nii');
        csf_file = dir('wc3*.nii');
        cd([sub_dir(nsub).folder,filesep,sub_dir(nsub).name]);
        cd('func');
        func_file = dir('wau*.nii');
        if isempty(func_file)
            func_file = dir('wu*.nii');
        end
        
        % Define paths to images
        func_img_path = [func_file(1).folder, filesep, func_file(1).name];  % Functional fMRI image
        wm_mask_path = [wm_file(1).folder, filesep, wm_file(1).name];  % White matter mask image
        csf_mask_path = [csf_file(1).folder, filesep, csf_file(1).name];         % CSF mask image

        % Load NIfTI files
        func_vol = spm_vol(func_img_path);
        wm_mask_vol = spm_vol(wm_mask_path);
        csf_mask_vol = spm_vol(csf_mask_path);

        % Read the image data
        func_data = spm_read_vols(func_vol);
        wm_mask_data = spm_read_vols(wm_mask_vol);
        csf_mask_data = spm_read_vols(csf_mask_vol);

        % Threshold masks (>0.9)
        wm_mask_thresh = wm_mask_data > 0.5;
        csf_mask_thresh = csf_mask_data > 0.5;

        % Reshape functional data for easier manipulation
        [nx, ny, nz, nt] = size(func_data);
        func_reshaped = reshape(func_data, [], nt);

        % Get indices of voxels in WM and CSF masks
        wm_voxel_indices = find(wm_mask_thresh);
        csf_voxel_indices = find(csf_mask_thresh);

        % Extract WM and CSF time series
        wm_time_series = func_reshaped(wm_voxel_indices, :);
        csf_time_series = func_reshaped(csf_voxel_indices, :);

        % Perform PCA on WM time series to extract first 5 PCs
        [~, wm_scores, ~] = pca(wm_time_series');
        wm_pcs = wm_scores(:, 1:5);  % First 5 principal components

        % Perform PCA on CSF time series to extract first 5 PCs
        [~, csf_scores, ~] = pca(csf_time_series');
        csf_pcs = csf_scores(:, 1:5);  % First 5 principal components

        % Display results
        fprintf('Extracted 5 principal components for White Matter and CSF signals.\n');
        
        % read movement regressor
        rp_file = dir('rp*.txt');
        movement_regressor = readmatrix([rp_file(1).folder,filesep,rp_file(1).name]);
        fprintf('Movement regerssor extracted.\n');
        
        %%% read ROI timeseries
        roi_time_series = zeros(718,nt);
        for ii = 1:length(roi_cortical_labels)
            label = roi_cortical_labels(ii);
            roi_mask = label_cortex == label;
            voxel_indices = find(roi_mask);

            roi_time_series(ii,:) = mean(func_reshaped(voxel_indices, :));
        end

        for jj = 1:length(roi_subcortical_labels_original)
            if any(roi_subcortical_labels_original(jj) == roi_subcortical_labels)
                label = roi_subcortical_labels_original(jj);
                roi_mask = label_subcortical == label;
                voxel_indices = find(roi_mask);

                roi_time_series(jj+360,:) = mean(func_reshaped(voxel_indices, :));
            end
        end
        fprintf('ROI timeseries extracted.\n');

        % Optional: Save WM and CSF PCs to a MAT file
        temp_timeseries_extracted{nsub} = {roi_time_series,movement_regressor',wm_pcs',csf_pcs'};
        temp_sub_list{nsub} = sub_dir(nsub).name;
    end
    AD_timeseries_extracted{ndir} = temp_timeseries_extracted;
    sub_list{ndir} = temp_sub_list;
end
%% 
cd(current_path);

save results/AD_timeseires_cortical_subcortical AD_timeseries_extracted sub_list