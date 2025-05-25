clear; clc; close all;
%% Load data
% Load Resting State Data
load('DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_subExclude_indiv_10_B.mat');
D_rest = D(2:end, tau ~= 0);
B_rest = B_mean(1:end, tau ~= 0);
sub_ids_rest = sub_ids(tau ~= 0);

% List of task files and names
task_names = {'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};
num_tasks = length(task_names);
D_tasks = cell(1, num_tasks);
B_tasks = cell(1, num_tasks);
sub_ids_tasks = cell(1, num_tasks);

% Load Task Data
for i = 1:num_tasks
    load(['DMs\DM_tfMRI_', task_names{i} ,'_cortical_subcortical_ext_fbDMD_noROInorm_subExclude_indiv_10_B.mat']);
    D_tasks{i} = D(2:end, tau ~= 0);
    B_tasks{i} = B(1:end, tau ~= 0);
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
B_rest = B_rest(:, idx_rest);
sub_ids_rest = sub_ids_rest(idx_rest);

% Subset Task Data to Common Subjects
for i = 1:num_tasks
    [~, idx_task] = ismember(common_sub_ids, sub_ids_tasks{i});
    D_tasks{i} = D_tasks{i}(:, idx_task);
    B_tasks{i} = B_tasks{i}(:, idx_task);
    sub_ids_tasks{i} = sub_ids_tasks{i}(idx_task);
end

%% Cosine distance
% Compute magnitude and phase of Rest Data
mag_rest = abs(D_rest(1:2:end,:));
phase_rest = angle(D_rest(1:2:end,:));
B_rest = B_rest(1:2:end,:);

% Compute magnitude and phase of Task Data
mag_tasks = cell(1, num_tasks);
phase_tasks = cell(1, num_tasks);

for i = 1:num_tasks
    mag_tasks{i} = abs(D_tasks{i}(1:2:end,:));
    phase_tasks{i} = angle(D_tasks{i}(1:2:end,:));
    B_tasks{i} = B_tasks{i}(1:2:end,:);
end

num_subjects = length(common_sub_ids);

% Initialize Cosine Similarity Matrices
% For Magnitude Comparisons
cosine_sim_mag_rest_tasks = zeros(num_subjects, num_tasks);
cosine_sim_mag_tasks = zeros(num_subjects, num_tasks, num_tasks);

% For Phase Comparisons
cosine_sim_phase_rest_tasks = zeros(num_subjects, num_tasks);
cosine_sim_phase_tasks = zeros(num_subjects, num_tasks, num_tasks);

% For Intensity Comparisons (B_tasks)
cosine_sim_B_rest_tasks = zeros(num_subjects, num_tasks);
cosine_sim_B_tasks = zeros(num_subjects, num_tasks, num_tasks);

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

% Compute Cosine Similarity Between Intensity (B_tasks) of Rest and Tasks
for i = 1:num_tasks
    X_B = B_rest;
    Y_B = B_tasks{i};
    numerator = sum(X_B .* Y_B, 1);
    denominator = sqrt(sum(X_B.^2, 1)) .* sqrt(sum(Y_B.^2, 1));
    cosine_similarity_B = numerator ./ denominator;
    cosine_sim_B_rest_tasks(:, i) = cosine_similarity_B;
end

% Compute Cosine Similarity Between Intensity (B_tasks) of Tasks
for i = 1:num_tasks
    for j = 1:num_tasks
        X_B = B_tasks{i};
        Y_B = B_tasks{j};
        numerator = sum(X_B .* Y_B, 1);
        denominator = sqrt(sum(X_B.^2, 1)) .* sqrt(sum(Y_B.^2, 1));
        cosine_similarity_B = numerator ./ denominator;
        cosine_sim_B_tasks(:, i, j) = cosine_similarity_B;
    end
end

% Mean Cosine Similarity for Magnitudes
mean_cosine_sim_mag_rest_tasks = mean(cosine_sim_mag_rest_tasks, 1);
ci_95_cosine_sim_mag_rest_tasks = 1.96 * std(cosine_sim_mag_rest_tasks) / sqrt(num_subjects);
mean_cosine_sim_mag_tasks = squeeze(mean(cosine_sim_mag_tasks, 1));

% Mean Cosine Similarity for Phases
mean_cosine_sim_phase_rest_tasks = mean(cosine_sim_phase_rest_tasks, 1);
ci_95_cosine_sim_phase_rest_tasks = 1.96 * std(cosine_sim_phase_rest_tasks) / sqrt(num_subjects);
mean_cosine_sim_phase_tasks = squeeze(mean(cosine_sim_phase_tasks, 1));

% Mean Cosine Similarity for Intensity (B_tasks)
mean_cosine_sim_B_rest_tasks = mean(cosine_sim_B_rest_tasks, 1);
ci_95_cosine_sim_B_rest_tasks = 1.96 * std(cosine_sim_B_rest_tasks) / sqrt(num_subjects);
mean_cosine_sim_B_tasks = squeeze(mean(cosine_sim_B_tasks, 1));

% Updated Task Names for Display
task_names_display = {'Emotional','Gambling','Language','Motor','Relational','Social','N-back'};

% Display Mean Cosine Similarity Between Intensity (B_tasks)
fprintf('Mean Cosine Similarity between Mode-Engagement Level of Rest and Tasks:\n');
for i = 1:num_tasks
    fprintf('Rest vs %s: %f\n', task_names_display{i}, mean_cosine_sim_B_rest_tasks(i));
end

fprintf('\nMean Cosine Similarity between Mode-Engagement Level of Tasks:\n');
B_similarity_table = array2table(mean_cosine_sim_B_tasks, 'VariableNames', task_names_display, 'RowNames', task_names_display);
disp(B_similarity_table);

% Display Mean Cosine Similarity Between Magnitudes
fprintf('Mean Cosine Similarity between Persistence Rates of Rest and Tasks:\n');
for i = 1:num_tasks
    fprintf('Rest vs %s: %f\n', task_names_display{i}, mean_cosine_sim_mag_rest_tasks(i));
end

fprintf('\nMean Cosine Similarity between Persistence Rates of Tasks:\n');
mag_similarity_table = array2table(mean_cosine_sim_mag_tasks, 'VariableNames', task_names_display, 'RowNames', task_names_display);
disp(mag_similarity_table);

% Display Mean Cosine Similarity Between Phases
fprintf('Mean Cosine Similarity between Progression Rates of Rest and Tasks:\n');
for i = 1:num_tasks
    fprintf('Rest vs %s: %f\n', task_names_display{i}, mean_cosine_sim_phase_rest_tasks(i));
end

fprintf('\nMean Cosine Similarity between Progression Rates of Tasks:\n');
phase_similarity_table = array2table(mean_cosine_sim_phase_tasks, 'VariableNames', task_names_display, 'RowNames', task_names_display);
disp(phase_similarity_table);

%% Cosine Distance Between Different Subjects in Same Task

% Initialize storage for Mean Cosine Similarities
mean_cosine_sim_mag_diff_subjects = zeros(1, num_tasks);
mean_cosine_sim_phase_diff_subjects = zeros(1, num_tasks);
mean_cosine_sim_B_diff_subjects = zeros(1, num_tasks); % New for B_tasks

% Initialize cell arrays to store similarity matrices
cosine_sim_mag_subjects = cell(1, num_tasks); % To store cosine similarity matrices for magnitudes
cosine_sim_phase_subjects = cell(1, num_tasks); % To store cosine similarity matrices for phases
cosine_sim_B_subjects = cell(1, num_tasks); % New for B_tasks

for i = 1:num_tasks
    %% Magnitude Data
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
    
    %% Phase Data
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
    
    %% Intensity Data (B_tasks) - New Section
    X_B = B_tasks{i}; % Size: [n_modes, n_subjects]
    % Normalize columns
    X_B_norm = sqrt(sum(X_B.^2,1));
    X_B_normalized = bsxfun(@rdivide, X_B, X_B_norm);
    % Compute cosine similarity matrix between subjects
    cosine_similarity_matrix_B = X_B_normalized' * X_B_normalized; % Size: [num_subjects, num_subjects]
    % Set diagonal to NaN to exclude self-similarity
    cosine_similarity_matrix_B(logical(eye(num_subjects))) = NaN;
    % Store the similarity matrix
    cosine_sim_B_subjects{i} = cosine_similarity_matrix_B;
    
    % Compute mean cosine similarity between different subjects
    mean_cosine_sim_B_diff_subjects(i) = nanmean(cosine_similarity_matrix_B(:));
end

% Display Mean Cosine Similarity Between Different Subjects for Each Task
fprintf('\nMean Cosine Similarity between Magnitudes of Different Subjects in Same Task:\n');
for i = 1:num_tasks
    fprintf('%s: %f\n', task_names_display{i}, mean_cosine_sim_mag_diff_subjects(i));
end

fprintf('\nMean Cosine Similarity between Phases of Different Subjects in Same Task:\n');
for i = 1:num_tasks
    fprintf('%s: %f\n', task_names_display{i}, mean_cosine_sim_phase_diff_subjects(i));
end

% Display Mean Cosine Similarity Between Intensity (B_tasks) of Different Subjects
fprintf('\nMean Cosine Similarity between Mode-Engagement Level of Different Subjects in Same Task:\n');
for i = 1:num_tasks
    fprintf('%s: %f\n', task_names_display{i}, mean_cosine_sim_B_diff_subjects(i));
end

%% Statistical tests
% Existing Statistical Tests for Magnitude and Phase
cosine_sim_mag_subjects_accum = [];
for i = 1:num_tasks
    cosine_similarity_matrix_mag = cosine_sim_mag_subjects{i};
    cosine_sim_mag_subjects_accum = [cosine_sim_mag_subjects_accum; cosine_similarity_matrix_mag(tril(ones(num_subjects,'logical'),-1))]; %#ok<AGROW>
end
cosine_sim_mag_tasks_accum = cosine_sim_mag_tasks(:, tril(ones(num_tasks,'logical'),-1));
[~,p_mag, ci_mag, stats_mag] = ttest2(cosine_sim_mag_subjects_accum, cosine_sim_mag_tasks_accum(:));
cohen_d_mag = (mean(cosine_sim_mag_subjects_accum) - mean(cosine_sim_mag_tasks_accum(:)))/stats_mag.sd;

cosine_sim_phase_subjects_accum = [];
for i = 1:num_tasks
    cosine_similarity_matrix_phase = cosine_sim_phase_subjects{i};
    cosine_sim_phase_subjects_accum = [cosine_sim_phase_subjects_accum; cosine_similarity_matrix_phase(tril(ones(num_subjects,'logical'),-1))]; %#ok<AGROW>
end
cosine_sim_phase_tasks_accum = cosine_sim_phase_tasks(:, tril(ones(num_tasks,'logical'),-1));
[~,p_phase, ci_phase, stats_phase] = ttest2(cosine_sim_phase_subjects_accum, cosine_sim_phase_tasks_accum(:));
cohen_d_phase = (mean(cosine_sim_phase_subjects_accum) - mean(cosine_sim_phase_tasks_accum(:)))/stats_phase.sd;

% New Statistical Tests for Intensity (B_tasks)
cosine_sim_B_subjects_accum = [];
for i = 1:num_tasks
    cosine_similarity_matrix_B = cosine_sim_B_subjects{i};
    cosine_sim_B_subjects_accum = [cosine_sim_B_subjects_accum; cosine_similarity_matrix_B(tril(ones(num_subjects,'logical'),-1))]; %#ok<AGROW>
end
cosine_sim_B_tasks_accum = cosine_sim_B_tasks(:, tril(ones(num_tasks,'logical'),-1));
[~,p_B, ci_B, stats_B] = ttest2(cosine_sim_B_subjects_accum, cosine_sim_B_tasks_accum(:));
cohen_d_B = (mean(cosine_sim_B_subjects_accum) - mean(cosine_sim_B_tasks_accum(:)))/stats_B.sd;

%% Compare Similarities
% Average similarity between different subjects in the same task
mean_mag_similarity_diff_subjects = mean(mean_cosine_sim_mag_diff_subjects);
mean_phase_similarity_diff_subjects = mean(mean_cosine_sim_phase_diff_subjects);
mean_B_similarity_diff_subjects = mean(mean_cosine_sim_B_diff_subjects); % New for B_tasks

% Average similarity between different tasks for the same subject
% For magnitudes
mean_mag_similarity_diff_tasks = mean(cosine_sim_mag_tasks(:, ~eye(num_tasks)), 'all'); % Exclude diagonal
% For phases
mean_phase_similarity_diff_tasks = mean(cosine_sim_phase_tasks(:, ~eye(num_tasks)), 'all');
% For Intensity (B_tasks) - New
mean_B_similarity_diff_tasks = mean(cosine_sim_B_tasks(:, ~eye(num_tasks)), 'all');

%%% Display Results
fprintf('\nAverage Cosine Similarity between Mode-Engagement Level of Different Subjects in Same Task:\n');
fprintf('Mean across all tasks: %f\n', mean_B_similarity_diff_subjects);

fprintf('\nAverage Cosine Similarity between Mode-Engagement Level of Different Tasks for Same Subject:\n');
fprintf('Mean across all task pairs: %f\n', mean_B_similarity_diff_tasks);
fprintf('cohen''s d : %f\n', cohen_d_B);
fprintf('t-value: %f\n', stats_B.tstat);
fprintf('p-value: %f\n', p_B);

fprintf('\nAverage Cosine Similarity between Persistence rate of Different Subjects in Same Task:\n');
fprintf('Mean across all tasks: %f\n', mean_mag_similarity_diff_subjects);

fprintf('\nAverage Cosine Similarity between Persistence rate of Different Tasks for Same Subject:\n');
fprintf('Mean across all task pairs: %f\n', mean_mag_similarity_diff_tasks);
fprintf('cohen''s d : %f\n', cohen_d_mag);
fprintf('t-value: %f\n', stats_mag.tstat);
fprintf('p-value: %f\n', p_mag);

fprintf('\nAverage Cosine Similarity between Progression rate of Different Subjects in Same Task:\n');
fprintf('Mean across all tasks: %f\n', mean_phase_similarity_diff_subjects);

fprintf('\nAverage Cosine Similarity between Progression rate of Different Tasks for Same Subject:\n');
fprintf('Mean across all task pairs: %f\n', mean_phase_similarity_diff_tasks);
fprintf('cohen''s d : %f\n', cohen_d_phase);
fprintf('t-value: %f\n', stats_phase.tstat);
fprintf('p-value: %f\n', p_phase);

