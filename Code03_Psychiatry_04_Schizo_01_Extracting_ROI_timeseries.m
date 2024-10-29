clear; clc;
current_path = pwd;

data_path = 'G:\COBRE_preprocessed\4197885';

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

%% functional image resampling
cd(data_path);
data_dir = dir('fmri*.nii');

reference_image_path = [current_path, filesep, 'atlas\MMP_cortex_resample_AD.nii'];

for nsub = 1:length(data_dir)
    cd([data_dir(nsub).folder,filesep,data_dir(nsub).name]);
    image_info = niftiinfo('file_s.nii');
    image_list = cell(image_info.ImageSize(4),1);
    for n_img = 1:image_info.ImageSize(4)
        image_list{n_img} = [data_dir(nsub).folder,filesep,data_dir(nsub).name,filesep,'file_s.nii,',num2str(n_img)];
    end
    
    matlabbatch = cell(1);
    matlabbatch{1}.spm.spatial.coreg.write.ref = {reference_image_path};    % Reference image
    matlabbatch{1}.spm.spatial.coreg.write.source = image_list;      % Source image to resample
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;        % Interpolation method (4th-degree B-spline)
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];    % No wrapping in any dimension
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;          % Do not mask the image
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';      % Prefix for the resampled image
    spm_jobman('run',matlabbatch);
    
end

%% ROI extraction
cd(data_path);
confound_dir = dir('fmri*.tsv');

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

schizo_timeseries_extracted = cell(length(data_dir),1);
schizo_sub_list = cell(length(data_dir),1);
for nsub = 1:length(data_dir)
    fprintf(['sub >> ',data_dir(nsub).name,'.\n']);
    cd([data_dir(nsub).folder,filesep,data_dir(nsub).name]);
    func_file = dir('rfile*.nii');
    if isempty(func_file)
        func_file = dir('wu*.nii');
    end
    
    func_img_path = [func_file(1).folder, filesep, func_file(1).name];  % Functional fMRI image

    % Load NIfTI files
    func_vol = spm_vol(func_img_path);

    % Read the image data
    func_data = spm_read_vols(func_vol);

    % Reshape functional data for easier manipulation
    [nx, ny, nz, nt] = size(func_data);
    func_reshaped = reshape(func_data, [], nt);
    
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
    
    cd([data_dir(nsub).folder,filesep,confound_dir(nsub).name]);
    
    confound_file = dir('fmri_*.tsv');
    
    confound_table = readtable(confound_file(1).name,'FileType', 'text', 'Delimiter', '\t');
    movement_regressor = confound_table(:,{'motion_tx','motion_ty','motion_tz','motion_rx','motion_ry','motion_rz'}).Variables;
    wm_pcs = confound_table(:,'wm_avg').Variables;
    csf_pcs = confound_table(:,'vent_avg').Variables;
    
    schizo_timeseries_extracted{nsub} = {roi_time_series,movement_regressor',wm_pcs',csf_pcs'};
    schizo_sub_list{nsub} = data_dir(nsub).name;
end
%% 
cd(current_path);

save results/Schizo_timeseires_cortical_subcortical schizo_timeseries_extracted schizo_sub_list