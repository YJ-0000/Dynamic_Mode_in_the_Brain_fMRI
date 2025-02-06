clear; clc;

%% Resting
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10_B
load secure_data/path_info;
gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
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

age = gene_data_table.Age_in_Yrs;
sex = strcmp(behav_data_table.Gender,'F');

age(tau == 0) = [];
sex(tau == 0) = [];

fprintf('=========== Resting ========= \n');
fprintf('Subject number: %d\n', length(age));
fprintf('Female: %d\n', sum(sex));
fprintf('age: (mean) %.2f  (std) %.2f \n', mean(age),std(age));


%% Tasks
fprintf('=========== Tasks ========= \n');
% List of task files and names
task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};
num_tasks = length(task_names);

% Load Task Data
for i = 1:num_tasks
    load(['DMs\DM_tfMRI_', task_names{i} ,'_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10.mat']);
    gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
    behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
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

    age = gene_data_table.Age_in_Yrs;
    sex = strcmp(behav_data_table.Gender,'F');

    age(tau == 0) = [];
    sex(tau == 0) = [];
    fprintf('-----> %s \n', task_names{i});
    fprintf('Subject number: %d\n', length(age));
    fprintf('Female: %d\n', sum(sex));
    fprintf('age: (mean) %.2f  (std) %.2f \n', mean(age),std(age));
    
end

%% Common subjects
% Load Resting State Data
load('DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10_B.mat');
sub_ids_rest = sub_ids(tau ~= 0);

% List of task files and names
task_names = {'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};
num_tasks = length(task_names);
sub_ids_tasks = cell(1, num_tasks);

% Load Task Data
for i = 1:num_tasks
    load(['DMs\DM_tfMRI_', task_names{i} ,'_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10_B.mat']);
    sub_ids_tasks{i} = sub_ids(tau ~= 0);
end

% Find Common Subjects Across All Datasets
common_sub_ids = sub_ids_rest;
for i = 1:num_tasks
    common_sub_ids = intersect(common_sub_ids, sub_ids_tasks{i});
end

% load age sex info
gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
for nrow = size(gene_data_table,1):-1:1
    if ~any(common_sub_ids==gene_data_table(nrow,'Subject').Variables)
        gene_data_table(nrow,:) = [];
    end
end
for nrow = size(behav_data_table,1):-1:1
    if ~any(common_sub_ids==behav_data_table(nrow,'Subject').Variables)
        behav_data_table(nrow,:) = [];
    end
end

gene_data_table = sortrows(gene_data_table, 'Subject');
behav_data_table = sortrows(behav_data_table, 'Subject');

age = gene_data_table.Age_in_Yrs;
sex = strcmp(behav_data_table.Gender,'F');

fprintf('=========== Common subjects (conducted all seven tasks) ========= \n');
fprintf('Subject number: %d\n', length(age));
fprintf('Female: %d\n', sum(sex));
fprintf('age: (mean) %.2f  (std) %.2f \n', mean(age),std(age));

