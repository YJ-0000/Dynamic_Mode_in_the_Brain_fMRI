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
%%% left hemisphere
% anterior group
cifti_data(cifti_data==361) = 361;
% Intralaminar group
cifti_data(cifti_data==362) = 362;
cifti_data(cifti_data==363) = 362;
cifti_data(cifti_data==364) = 362;
cifti_data(cifti_data==373) = 362;
% Lateral group
cifti_data(cifti_data==365) = 363;
cifti_data(cifti_data==367) = 363;
% Posterior group
cifti_data(cifti_data==366) = 364;
cifti_data(cifti_data==368) = 364;
cifti_data(cifti_data==371) = 364;
cifti_data(cifti_data==374) = 364;
cifti_data(cifti_data==375) = 364;
cifti_data(cifti_data==376) = 364;
cifti_data(cifti_data==377) = 364;
% Medial group
cifti_data(cifti_data==369) = 365;
cifti_data(cifti_data==370) = 365;
cifti_data(cifti_data==372) = 365;
% Ventral group
cifti_data(cifti_data==378) = 366;
cifti_data(cifti_data==379) = 366;
cifti_data(cifti_data==380) = 366;
cifti_data(cifti_data==381) = 366;

% other subcortical
for ii = 1:12
    cifti_data(cifti_data==(381+ii)) = 366+ii;
end


%%% right hemisphere
% anterior group
cifti_data(cifti_data==361+33) = 361+33;
% Intralaminar group
cifti_data(cifti_data==362+33) = 362+33;
cifti_data(cifti_data==363+33) = 362+33;
cifti_data(cifti_data==364+33) = 362+33;
cifti_data(cifti_data==373+33) = 362+33;
% Lateral group
cifti_data(cifti_data==365+33) = 363+33;
cifti_data(cifti_data==367+33) = 363+33;
% Posterior group
cifti_data(cifti_data==366+33) = 364+33;
cifti_data(cifti_data==368+33) = 364+33;
cifti_data(cifti_data==371+33) = 364+33;
cifti_data(cifti_data==374+33) = 364+33;
cifti_data(cifti_data==375+33) = 364+33;
cifti_data(cifti_data==376+33) = 364+33;
cifti_data(cifti_data==377+33) = 364+33;
% Medial group
cifti_data(cifti_data==369+33) = 365+33;
cifti_data(cifti_data==370+33) = 365+33;
cifti_data(cifti_data==372+33) = 365+33;
% Ventral group
cifti_data(cifti_data==378+33) = 366+33;
cifti_data(cifti_data==379+33) = 366+33;
cifti_data(cifti_data==380+33) = 366+33;
cifti_data(cifti_data==381+33) = 366+33;


% other subcortical
for ii = 1:12
    cifti_data(cifti_data==(381+33+ii)) = 366+33+ii;
end

% right subcortical
% other subcortical
for ii = 1:18
    cifti_data(cifti_data==(393+ii)) = 378+ii;
end

% cerebellum
for ii = 1:24
    cifti_data(cifti_data==(426+ii)) = 396+ii;
end
%%
cifti_data(cifti_data<=360) = 0;
%%
for nroi = 1:18
    disp([360+nroi, sum(cifti_data==(360+nroi)), 360+18+nroi, sum(cifti_data==(360+18+nroi))]);
end
for nroi = 1:12
    disp([396+nroi, sum(cifti_data==(396+nroi)), 396+12+nroi, sum(cifti_data==(396+12+nroi))]);
end

%%
labels = cifti_struct_create_from_template(cifti, cifti_data, 'dscalar');
cifti_write(labels, 'atlas/Subcortical_6thalamus.dscalar.nii');