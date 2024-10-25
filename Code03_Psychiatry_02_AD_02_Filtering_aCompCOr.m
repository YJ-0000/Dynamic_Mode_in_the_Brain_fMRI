clear; clc;
current_path = pwd;

load('results/AD_timeseires_cortical_subcortical.mat');

data_path = 'D:\ADNI';

%%
cd(data_path);
data_dir = dir;
data_dir(~[data_dir.isdir]) = [];
data_dir(1:2) = [];
data_dir_name = {data_dir.name};

excel_file_path = [data_path,filesep,'ADNI_prep_list.xlsx'];

for i = 1:length(data_dir_name)
    sheet_name = data_dir_name{i};
    fprintf('Reading sheet: %s\n', sheet_name);
    data = readtable(excel_file_path, 'Sheet', sheet_name);
end

cd(current_path);
%%
% Define the frequency range for the band-pass filter
lowCutoff = 0.01;  % Hz
highCutoff = 0.1;  % Hz

%%
AD_timeseries_extracted_denoised_filtered = cell(size(AD_timeseries_extracted));
AD_TR_info = cell(size(AD_timeseries_extracted));
AD_sub_info = cell(size(AD_timeseries_extracted));
for ndir = 1:length(AD_timeseries_extracted)
    sub_info = readtable(excel_file_path, 'Sheet', data_dir_name{ndir});
    AD_sub_info{ndir} = sub_info;
    fprintf('Reading sheet: %s\n', sheet_name);
    temp_timeseries_extracted = AD_timeseries_extracted{ndir};
    temp_sub_list = sub_list{ndir};
    
    temp_timeseries_extracted_denoised_filtered = cell(size(temp_timeseries_extracted));
    temp_TR_info = zeros(size(temp_timeseries_extracted));
    for nsub = 1:length(temp_timeseries_extracted)
        sub_timeseries_extracted = temp_timeseries_extracted{nsub};
        sub_name = temp_sub_list{nsub};
        fprintf(['sub >> ',sub_name,'.\n']);
        
        sub_idx = find(strcmp(sub_info.subject_id',sub_name));
        assert(~isempty(sub_idx));
        
        tr = sub_info(sub_idx,'func_TR').Variables;
        temp_TR_info(nsub) = tr;
        
        Fs = 1/tr;
        
        fmri_data = sub_timeseries_extracted{1}';
        
        len_time = size(fmri_data,1);
        
        move_regressor = sub_timeseries_extracted{2}';
        move_regressor = [move_regressor, [zeros(1,6);(move_regressor(2:end,:)-move_regressor(1:end-1,:))]];
        move_regressor = [move_regressor, move_regressor.^2];
        
        noise_mat = [ones(len_time,1), (1:len_time)',((1:len_time).^2)',move_regressor, sub_timeseries_extracted{3}', sub_timeseries_extracted{4}'];
        
        fmri_data2 = (fmri_data - noise_mat*pinv(noise_mat)*(fmri_data));
        fmri_data2 = bandpass(fmri_data2, [lowCutoff,highCutoff], Fs)';
        
        temp_timeseries_extracted_denoised_filtered{nsub} = fmri_data2;
        
        fprintf('Filtered and aCompCor.\n');
    end
    AD_timeseries_extracted_denoised_filtered{ndir} = temp_timeseries_extracted_denoised_filtered;
    AD_TR_info{ndir} = temp_TR_info;
end

%%
save results/AD_timeseires_cortical_subcortical_denoised_filtered sub_list AD_timeseries_extracted_denoised_filtered AD_sub_info