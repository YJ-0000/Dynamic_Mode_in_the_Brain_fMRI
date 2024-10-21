clear; clc; close all;
%% Load data
% Load Resting State Data
load('DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10.mat');
D_rest = D(2:end, tau ~= 0).^(0.72/1.5);
sub_ids_rest = sub_ids(tau ~= 0);

% List of task files and names
task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};
num_tasks = length(task_names);
D_tasks = cell(1, num_tasks);
sub_ids_tasks = cell(1, num_tasks);

% Load Task Data
for i = 1:num_tasks
    load(['DMs\DM_tfMRI_', task_names{i} ,'_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10.mat']);
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

%% Cosine distance
% Compute magnitude and phase of Rest Data
mag_rest = abs(D_rest(1:2:end,:));
phase_rest = angle(D_rest(1:2:end,:));

% Compute magnitude and phase of Task Data
mag_tasks = cell(1, num_tasks);
phase_tasks = cell(1, num_tasks);

for i = 1:num_tasks
    mag_tasks{i} = abs(D_tasks{i}(1:2:end,:));
    phase_tasks{i} = angle(D_tasks{i}(1:2:end,:));
end

num_subjects = length(common_sub_ids);

% For Magnitude Comparisons
cosine_sim_mag_rest_tasks = zeros(num_subjects, num_tasks);
cosine_sim_mag_tasks = zeros(num_subjects, num_tasks, num_tasks);

% For Phase Comparisons
cosine_sim_phase_rest_tasks = zeros(num_subjects, num_tasks);
cosine_sim_phase_tasks = zeros(num_subjects, num_tasks, num_tasks);

% Compute Cosine Similarity Between Magnitudes of Rest and Tasks
for i = 1:num_tasks
    X_mag = mag_rest;
    Y_mag = mag_tasks{i};
    numerator = sum(X_mag .* Y_mag, 1);
    denominator = sqrt(sum(X_mag.^2, 1)) .* sqrt(sum(Y_mag.^2, 1));
    cosine_similarity_mag = numerator ./ denominator;
    cosine_sim_mag_rest_tasks(:, i) = cosine_similarity_mag;
end

% Compute Cosine Similarity Between Magnitudes of Tasks
for i = 1:num_tasks
    for j = 1:num_tasks
        X_mag = mag_tasks{i};
        Y_mag = mag_tasks{j};
        numerator = sum(X_mag .* Y_mag, 1);
        denominator = sqrt(sum(X_mag.^2, 1)) .* sqrt(sum(Y_mag.^2, 1));
        cosine_similarity_mag = numerator ./ denominator;
        cosine_sim_mag_tasks(:, i, j) = cosine_similarity_mag;
    end
end

% Compute Cosine Similarity Between Phases of Rest and Tasks
for i = 1:num_tasks
    X_phase = phase_rest;
    Y_phase = phase_tasks{i};
    numerator = sum(X_phase .* Y_phase, 1);
    denominator = sqrt(sum(X_phase.^2, 1)) .* sqrt(sum(Y_phase.^2, 1));
    cosine_similarity_phase = numerator ./ denominator;
    cosine_sim_phase_rest_tasks(:, i) = cosine_similarity_phase;
end

% Compute Cosine Similarity Between Phases of Tasks
for i = 1:num_tasks
    for j = 1:num_tasks
        X_phase = phase_tasks{i};
        Y_phase = phase_tasks{j};
        numerator = sum(X_phase .* Y_phase, 1);
        denominator = sqrt(sum(X_phase.^2, 1)) .* sqrt(sum(Y_phase.^2, 1));
        cosine_similarity_phase = numerator ./ denominator;
        cosine_sim_phase_tasks(:, i, j) = cosine_similarity_phase;
    end
end

% Mean Cosine Similarity for Magnitudes
mean_cosine_sim_mag_rest_tasks = mean(cosine_sim_mag_rest_tasks, 1);
mean_cosine_sim_mag_tasks = squeeze(mean(cosine_sim_mag_tasks, 1));

% Mean Cosine Similarity for Phases
mean_cosine_sim_phase_rest_tasks = mean(cosine_sim_phase_rest_tasks, 1);
mean_cosine_sim_phase_tasks = squeeze(mean(cosine_sim_phase_tasks, 1));

% Display Mean Cosine Similarity Between Magnitudes
fprintf('Mean Cosine Similarity between Magnitudes of Rest and Tasks:\n');
for i = 1:num_tasks
    fprintf('Rest vs %s: %f\n', task_names{i}, mean_cosine_sim_mag_rest_tasks(i));
end

fprintf('\nMean Cosine Similarity between Magnitudes of Tasks:\n');
mag_similarity_table = array2table(mean_cosine_sim_mag_tasks, 'VariableNames', task_names, 'RowNames', task_names);
disp(mag_similarity_table);

% Display Mean Cosine Similarity Between Phases
fprintf('Mean Cosine Similarity between Phases of Rest and Tasks:\n');
for i = 1:num_tasks
    fprintf('Rest vs %s: %f\n', task_names{i}, mean_cosine_sim_phase_rest_tasks(i));
end

fprintf('\nMean Cosine Similarity between Magnitudes of Tasks:\n');
phase_similarity_table = array2table(mean_cosine_sim_phase_tasks, 'VariableNames', task_names, 'RowNames', task_names);
disp(phase_similarity_table);


%% Cosine Distance Between Different Subjects in Same Task

% For each task, compute cosine similarities between different subjects
cosine_sim_mag_subjects = cell(1, num_tasks); % To store cosine similarity matrices for magnitudes
cosine_sim_phase_subjects = cell(1, num_tasks); % To store cosine similarity matrices for phases

for i = 1:num_tasks
    % Magnitude data
    X_mag = mag_tasks{i}; % Size: [n_modes, n_subjects]
    % Normalize columns
    X_mag_norm = sqrt(sum(X_mag.^2,1));
    X_mag_normalized = bsxfun(@rdivide, X_mag, X_mag_norm);
    % Compute cosine similarity matrix between subjects
    cosine_similarity_matrix_mag = X_mag_normalized' * X_mag_normalized; % Size: [num_subjects, num_subjects]
    % Set diagonal to NaN to exclude self-similarity
    cosine_similarity_matrix_mag(logical(eye(num_subjects))) = NaN;
    % Store the similarity matrix
    cosine_sim_mag_subjects{i} = cosine_similarity_matrix_mag;
    
    % Compute mean cosine similarity between different subjects
    mean_cosine_sim_mag_diff_subjects(i) = nanmean(cosine_similarity_matrix_mag(:));
    
    % Phase data
    X_phase = phase_tasks{i}; % Size: [n_modes, n_subjects]
    % Normalize columns
    X_phase_norm = sqrt(sum(X_phase.^2,1));
    X_phase_normalized = bsxfun(@rdivide, X_phase, X_phase_norm);
    % Compute cosine similarity matrix between subjects
    cosine_similarity_matrix_phase = X_phase_normalized' * X_phase_normalized; % Size: [num_subjects, num_subjects]
    % Set diagonal to NaN to exclude self-similarity
    cosine_similarity_matrix_phase(logical(eye(num_subjects))) = NaN;
    % Store the similarity matrix
    cosine_sim_phase_subjects{i} = cosine_similarity_matrix_phase;
    
    % Compute mean cosine similarity between different subjects
    mean_cosine_sim_phase_diff_subjects(i) = nanmean(cosine_similarity_matrix_phase(:));
end

% Display Mean Cosine Similarity Between Different Subjects for Each Task
fprintf('\nMean Cosine Similarity between Magnitudes of Different Subjects in Same Task:\n');
for i = 1:num_tasks
    fprintf('%s: %f\n', task_names{i}, mean_cosine_sim_mag_diff_subjects(i));
end

fprintf('\nMean Cosine Similarity between Phases of Different Subjects in Same Task:\n');
for i = 1:num_tasks
    fprintf('%s: %f\n', task_names{i}, mean_cosine_sim_phase_diff_subjects(i));
end

%% Statistical tests
cosine_sim_mag_subjects_accum = [];
for i = 1:num_tasks
    cosine_similarity_matrix_mag = cosine_sim_mag_subjects{i};
    cosine_sim_mag_subjects_accum = [cosine_sim_mag_subjects_accum;cosine_similarity_matrix_mag(~isnan(cosine_similarity_matrix_mag))];
end
cosine_sim_mag_tasks_accum = cosine_sim_mag_tasks(:,~eye(num_tasks));
[~,p_mag] = ttest2(cosine_sim_mag_subjects_accum,cosine_sim_mag_tasks_accum(:));

cosine_sim_phase_subjects_accum = [];
for i = 1:num_tasks
    cosine_similarity_matrix_phase = cosine_sim_phase_subjects{i};
    cosine_sim_phase_subjects_accum = [cosine_sim_phase_subjects_accum;cosine_similarity_matrix_phase(~isnan(cosine_similarity_matrix_phase))];
end
cosine_sim_phase_tasks_accum = cosine_sim_phase_tasks(:,~eye(num_tasks));
[~,p_phase] = ttest2(cosine_sim_phase_subjects_accum,cosine_sim_phase_tasks_accum(:));

%% Compare Similarities
% Average similarity between different subjects in the same task
mean_mag_similarity_diff_subjects = mean(mean_cosine_sim_mag_diff_subjects);
mean_phase_similarity_diff_subjects = mean(mean_cosine_sim_phase_diff_subjects);

% Average similarity between different tasks for the same subject
% For magnitudes
mean_mag_similarity_diff_tasks = mean(cosine_sim_mag_tasks(:,~eye(num_tasks)),'all'); % Exclude diagonal
% For phases
mean_phase_similarity_diff_tasks = mean(cosine_sim_phase_tasks(:,~eye(num_tasks)),'all');

%%% Display Results
fprintf('\nAverage Cosine Similarity between Magnitudes of Different Subjects in Same Task:\n');
fprintf('Mean across all tasks: %f\n', mean_mag_similarity_diff_subjects);

fprintf('\nAverage Cosine Similarity between Magnitudes of Different Tasks for Same Subject:\n');
fprintf('Mean across all task pairs: %f\n', mean_mag_similarity_diff_tasks);
fprintf('p-value: %f\n', p_mag);

fprintf('\nAverage Cosine Similarity between Phases of Different Subjects in Same Task:\n');
fprintf('Mean across all tasks: %f\n', mean_phase_similarity_diff_subjects);

fprintf('\nAverage Cosine Similarity between Phases of Different Tasks for Same Subject:\n');
fprintf('Mean across all task pairs: %f\n', mean_phase_similarity_diff_tasks);
fprintf('p-value: %f\n', p_phase);
