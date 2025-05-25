clear; clc; close all;

test_for = 'phase'; % magnitude or phase

%% Load data
% Load Resting State Data
load('DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_subExclude_indiv_10_B.mat');
D_rest = D(2:end, tau ~= 0);
sub_ids_rest = sub_ids(tau ~= 0);

% List of task files and names
task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};
num_tasks = length(task_names);
D_tasks = cell(1, num_tasks);
sub_ids_tasks = cell(1, num_tasks);

% Load Task Data
for i = 1:num_tasks
    load(['DMs\DM_tfMRI_', task_names{i} ,'_cortical_subcortical_ext_fbDMD_noROInorm_subExclude_indiv_10_B.mat']);
    D_tasks{i} = D(2:end, tau ~= 0);
    sub_ids_tasks{i} = sub_ids(tau ~= 0);
end

% Find Common Subjects Across All Datasets
common_sub_ids = sub_ids_rest;
for i = 1:num_tasks
    common_sub_ids = intersect(common_sub_ids, sub_ids_tasks{i});
end

% Subset Rest Data to Common Subjects
[~, idx_rest] = ismember(common_sub_ids, sub_ids_rest);
D_rest = D_rest(:, idx_rest);
sub_ids_rest = sub_ids_rest(idx_rest);

% Subset Task Data to Common Subjects
for i = 1:num_tasks
    [~, idx_task] = ismember(common_sub_ids, sub_ids_tasks{i});
    D_tasks{i} = D_tasks{i}(:, idx_task);
    sub_ids_tasks{i} = sub_ids_tasks{i}(idx_task);
end

%% t-test
idx_order = [1,3,9,7,5];

phase_rest = angle(D_rest(idx_order,:));
phase_tasks = cell(1, num_tasks);
for i = 1:num_tasks
    phase_tasks{i} = angle(D_tasks{i}(idx_order,:));
end

for n_dm = 1:5
    [h,p] = ttest(phase_rest(n_dm,:));
    fprintf('Significance of DM%d: mean theta = %.6f , p=%.6f \n',n_dm, mean(phase_rest(n_dm,:)), p);
end
disp(' ');

for n_task = 1:num_tasks
    phase_task_1 = phase_tasks{n_task};
    for n_dm = 1:5
        [h,p,~,stat] = ttest(phase_task_1(n_dm,:));
        cohen_d = mean(phase_task_1(n_dm,:)) / stat.sd;
        fprintf('Significance of DM%d during %s: cohen''s d=%.6f , t=%.6f,  p=%.11f \n', ...
            n_dm,task_names{n_task}, cohen_d, stat.tstat, p);
    end
    disp(' ');
end