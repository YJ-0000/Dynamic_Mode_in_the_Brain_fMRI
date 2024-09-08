clear; clc;

%%
% Initialize SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%%
% Specify the path to your NIfTI file
nifti_file = 'atlas/HCPex.nii';

% Load the NIfTI header
V = spm_vol(nifti_file);
[atlas_data_HCPex, XYZ_HCPex] = spm_read_vols(V);
labels_HCPex = unique(atlas_data_HCPex(:));
labels_HCPex(labels_HCPex == 0) = [];

% Specify the path to your NIfTI file
nifti_file = 'atlas/Cerebellum_resample.nii';

% Load the NIfTI header
V = spm_vol(nifti_file);
[atlas_data_Cerebellum, XYZ_Cerebellum] = spm_read_vols(V);
atlas_data_Cerebellum = round(atlas_data_Cerebellum);
labels_Cerebellum = unique(atlas_data_Cerebellum(:));
labels_Cerebellum(labels_Cerebellum == 0) = [];
for ii = 1:length(labels_Cerebellum)
    atlas_data_Cerebellum(atlas_data_Cerebellum==labels_Cerebellum(ii)) = max(labels_HCPex)+ii;
    labels_Cerebellum(ii) = max(labels_HCPex)+ii;
end

atlas_data_all = atlas_data_HCPex;
atlas_data_all(atlas_data_all == 0) = atlas_data_Cerebellum(atlas_data_all == 0);
XYZ = XYZ_HCPex;

%%
cifti = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
cifti_data = zeros(size(cifti.cdata));
for n_model = 3:length(cifti.diminfo{1, 1}.models)

    disp(cifti.diminfo{1, 1}.models{1, n_model}.count);
    
    cifti_MNI_coords = cifti.diminfo{1, 1}.vol.sform  * [cifti.diminfo{1, 1}.models{1, n_model}.voxlist;ones(1,cifti.diminfo{1, 1}.models{1, n_model}.count)];
    cifti_MNI_coords = cifti_MNI_coords(1:3,:);
    
    temp_data = zeros(cifti.diminfo{1, 1}.models{1, n_model}.count,1);
    parfor n_voxel = 1:cifti.diminfo{1, 1}.models{1, n_model}.count
        current_cifti_MNI_coords = cifti_MNI_coords(:,n_voxel);
        
        dist_matrix = sum((XYZ - current_cifti_MNI_coords).^2,1);
        [min_dixt,min_voxel_n] = min(dist_matrix);
        
        temp_data(n_voxel) = atlas_data_all(min_voxel_n);
    end
    cifti_data(cifti.diminfo{1, 1}.models{1, n_model}.start:cifti.diminfo{1, 1}.models{1, n_model}.start+cifti.diminfo{1, 1}.models{1, n_model}.count-1) = temp_data;
    disp(cifti.diminfo{1, 1}.models{1, n_model}.struct);
end

%%
labels = cifti_struct_create_from_template(cifti, cifti_data, 'dscalar');
cifti_write(labels, 'atlas/Cortical_Subcortical.dscalar.nii');