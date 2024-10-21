clear; clc; close all;

test_for = 'phase'; % magnitude or phase

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
%% coefficients
idx_order = [1,5,9,7,3];

figure; 
subplot(2,4,1);
bar(mean(abs(D_rest(idx_order,:)),2));
title('REST');
ylim([0.94,1]);
for ii = 1:7
    subplot(2,4,ii+1);
    bar(mean(abs(D_tasks{ii}(idx_order,:)),2));
    title(task_names{ii});
    ylim([0.94,1]);
end

%% t-test
% Normalize D
% normalized_D_rest = abs(D_rest(1:2:end,:)) ./ mean(abs(D_rest(1:2:end,:)), 2);
% normalized_D_tasks = cell(1, num_tasks);
% for i = 1:num_tasks
%     normalized_D_tasks{i} = abs(D_tasks{i}(1:2:end,:)) ./ mean(abs(D_tasks{i}(1:2:end,:)), 2);
% end

if strcmp(test_for,'magnitude')
    normalized_D_rest = abs(D_rest(idx_order,:)) ;
    normalized_D_tasks = cell(1, num_tasks);
    for i = 1:num_tasks
        normalized_D_tasks{i} = abs(D_tasks{i}(idx_order,:)) ;
    end
elseif strcmp(test_for,'phase')
    normalized_D_rest = angle(D_rest(idx_order,:));
    normalized_D_tasks = cell(1, num_tasks);
    for i = 1:num_tasks
        normalized_D_tasks{i} = angle(D_tasks{i}(idx_order,:));
    end
else
    error('Magnitude or phase??')
end

% Number of modes
n_modes = size(normalized_D_rest, 1);

% Initialize variables to store t-test results
p_values_rest_vs_tasks = zeros(n_modes, num_tasks);
t_stats_rest_vs_tasks = zeros(n_modes, num_tasks);

p_values_tasks_vs_tasks = zeros(n_modes, num_tasks, num_tasks);
t_stats_tasks_vs_tasks = zeros(n_modes, num_tasks, num_tasks);

% Perform t-test between Rest and each Task for each mode
for mode_idx = 1:n_modes
    rest_data = normalized_D_rest(mode_idx, :);
    for task_idx = 1:num_tasks
        task_data = normalized_D_tasks{task_idx}(mode_idx, :);
        [~, p, ~, stats] = ttest(rest_data, task_data);
        p_values_rest_vs_tasks(mode_idx, task_idx) = p;
        t_stats_rest_vs_tasks(mode_idx, task_idx) = stats.tstat;
    end
end

% Perform t-test between each pair of Tasks for each mode
for mode_idx = 1:n_modes
    for task_idx1 = 1:num_tasks
        data1 = normalized_D_tasks{task_idx1}(mode_idx, :);
        for task_idx2 = task_idx1+1:num_tasks
            data2 = normalized_D_tasks{task_idx2}(mode_idx, :);
            [~, p, ~, stats] = ttest(data1, data2);
            p_values_tasks_vs_tasks(mode_idx, task_idx1, task_idx2) = p;
            t_stats_tasks_vs_tasks(mode_idx, task_idx1, task_idx2) = stats.tstat;
        end
    end
end

% Correct for multiple comparisons using False Discovery Rate (FDR)
% Concatenate all p-values into a single vector for FDR correction
all_p_values = [p_values_rest_vs_tasks(:); p_values_tasks_vs_tasks(p_values_tasks_vs_tasks > 0)];

% Apply FDR correction
[~, ~, ~, adj_p_values] = fdr_bh(all_p_values);

% Split adjusted p-values back into their respective matrices
p_values_rest_vs_tasks_adj = reshape(adj_p_values(1:numel(p_values_rest_vs_tasks)), size(p_values_rest_vs_tasks));
p_values_tasks_vs_tasks_adj = zeros(size(p_values_tasks_vs_tasks));
idx = numel(p_values_rest_vs_tasks) + 1;
for task_idx1 = 1:num_tasks
    for task_idx2 = task_idx1+1:num_tasks
        num_p = n_modes;
        p_values_tasks_vs_tasks_adj(:, task_idx1, task_idx2) = adj_p_values(idx:idx+num_p-1);
        idx = idx + num_p;
    end
end

% Set non-computed entries to NaN
p_values_tasks_vs_tasks_adj(p_values_tasks_vs_tasks_adj == 0) = NaN;

% Bonferroni correction
% Number of comparisons
N_rest_vs_tasks = n_modes * num_tasks;
N_tasks_vs_tasks = n_modes * nchoosek(num_tasks, 2);

% Bonferroni-adjusted p-values
p_values_rest_vs_tasks_bonf = p_values_rest_vs_tasks * N_rest_vs_tasks;
p_values_rest_vs_tasks_bonf(p_values_rest_vs_tasks_bonf > 1) = 1;

p_values_tasks_vs_tasks_bonf = p_values_tasks_vs_tasks * N_tasks_vs_tasks;
p_values_tasks_vs_tasks_bonf(p_values_tasks_vs_tasks_bonf > 1) = 1;
p_values_tasks_vs_tasks_bonf(p_values_tasks_vs_tasks_bonf == 0) = NaN;

%% Visualization using heatmap
alpha = 0.05; % Significance level

n = 256; % Number of colors
half = floor(n/2);
    
% Blue to white
cmap1 = [linspace(0,1,half)', linspace(0,1,half)', linspace(1,1,half)'];

% White to red
cmap2 = [linspace(1,1,n-half)', linspace(1,0,n-half)', linspace(1,0,n-half)'];

% Combine to create red-white-blue colormap
redWhiteBlue = [cmap1; cmap2];

% Rest vs Tasks Heatmap
figure;
% Create t-value heatmap
imagesc(-t_stats_rest_vs_tasks);
colormap(redWhiteBlue);
colorbar;
mode_names = {'Principal','CEN-to-DMN','DMN-to-CEN','Bi-asym','FV-to-S'};
caxis([-max(abs(t_stats_rest_vs_tasks(:))), max(abs(t_stats_rest_vs_tasks(:)))]);
set(gca, 'XTick', 1:num_tasks, 'XTickLabel', task_names, 'YTick', 1:n_modes,'YTickLabel', mode_names);
% xlabel('Tasks');
% ylabel('Modes');
title('T-values: Rest vs Tasks');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% Overlay significance markers
hold on;
for mode_idx = 1:n_modes
    for task_idx = 1:num_tasks
        x = task_idx;
        y = mode_idx;
        p_fdr = p_values_rest_vs_tasks_adj(mode_idx, task_idx);
        p_bonf = p_values_rest_vs_tasks_bonf(mode_idx, task_idx);
        if p_bonf < alpha
            text(x, y, '*', 'HorizontalAlignment', 'Center', 'Color', 'k', 'FontSize', 14);
        elseif p_fdr < alpha
            text(x, y, '+', 'HorizontalAlignment', 'Center', 'Color', 'k', 'FontSize', 14);
        end
    end
end
hold off;

% Tasks vs Tasks Heatmap
% Prepare data for heatmap
n_task_pairs = nchoosek(num_tasks, 2);
task_pairs = nchoosek(1:num_tasks, 2);

t_values_tasks_vs_tasks = zeros(n_modes, n_task_pairs);
p_values_tasks_vs_tasks_adj_flat = zeros(n_modes, n_task_pairs);
p_values_tasks_vs_tasks_bonf_flat = zeros(n_modes, n_task_pairs);
pair_labels = cell(1, n_task_pairs);

for pair_idx = 1:n_task_pairs
    idx1 = task_pairs(pair_idx, 1);
    idx2 = task_pairs(pair_idx, 2);
    t_values_tasks_vs_tasks(:, pair_idx) = t_stats_tasks_vs_tasks(:, idx1, idx2);
    p_values_tasks_vs_tasks_adj_flat(:, pair_idx) = p_values_tasks_vs_tasks_adj(:, idx1, idx2);
    p_values_tasks_vs_tasks_bonf_flat(:, pair_idx) = p_values_tasks_vs_tasks_bonf(:, idx1, idx2);
    pair_labels{pair_idx} = [task_names{idx1} ' vs ' task_names{idx2}];
end

% Tasks vs Tasks Heatmaps (Square Matrices)
% For each mode, create a heatmap of t-values between tasks
num_figures = 5; % Number of figures to display
modes_per_figure = ceil(n_modes / num_figures);

for fig_idx = 1:num_figures
    figure;
    mode_start = (fig_idx - 1) * modes_per_figure + 1;
    mode_end = min(fig_idx * modes_per_figure, n_modes);
    num_modes_in_fig = mode_end - mode_start + 1;
    
    for subplot_idx = 1:num_modes_in_fig
        mode_idx = mode_start + subplot_idx - 1;
        t_matrix = squeeze(t_stats_tasks_vs_tasks(mode_idx, :, :));
        p_fdr_matrix = squeeze(p_values_tasks_vs_tasks_adj(mode_idx, :, :));
        p_bonf_matrix = squeeze(p_values_tasks_vs_tasks_bonf(mode_idx, :, :));
        
        subplot(ceil(sqrt(num_modes_in_fig)), ceil(sqrt(num_modes_in_fig)), subplot_idx);
        imagesc(t_matrix);
        colormap(flipud(jet));
        caxis([-max(abs(t_matrix(:))), max(abs(t_matrix(:)))]);
        colormap(redWhiteBlue);
        colorbar;
        title(['Mode ' num2str(mode_idx)]);
        set(gca, 'XTick', 1:num_tasks, 'XTickLabel', task_names, 'XTickLabelRotation', 45);
        set(gca, 'YTick', 1:num_tasks, 'YTickLabel', task_names);
        
        % Overlay significance markers
        hold on;
        for task_idx1 = 1:num_tasks
            for task_idx2 = 1:num_tasks
                x = task_idx2;
                y = task_idx1;
                p_fdr = p_fdr_matrix(task_idx1, task_idx2);
                p_bonf = p_bonf_matrix(task_idx1, task_idx2);
                if p_bonf < alpha
                    text(x, y, '*', 'HorizontalAlignment', 'Center', 'Color', 'k', 'FontSize', 12);
                elseif p_fdr < alpha
                    text(x, y, '+', 'HorizontalAlignment', 'Center', 'Color', 'k', 'FontSize', 12);
                end
            end
        end
        hold off;
    end
    title(['Tasks vs Tasks (Modes ' num2str(mode_start) ' to ' num2str(mode_end) ')']);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
    
end

%% Function: FDR Correction
function [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report)
% Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR).
% This function adjusts p-values using the FDR method.

% Default values for optional parameters
if nargin<2
    q=.05;
end
if nargin<3
    method='pdep'; % Dependent or independent
end
if nargin<4
    report='no';
end

pvals=pvals(:);
m=length(pvals);

if strcmpi(method,'pdep')
    [sorted_pvals,sorted_index]=sort(pvals);
    [~,unsort_order]=sort(sorted_index);
    thresh=(1:m)'/m*q;
    wtd_p=sorted_pvals./thresh;
    max_id=find(sorted_pvals<=thresh, 1, 'last' );
elseif strcmpi(method,'dep')
    % This is the "two-stage" procedure
    [sorted_pvals,sorted_index]=sort(pvals);
    [~,unsort_order]=sort(sorted_index);
    c=sum(1./(1:m));
    thresh=(1:m)'/m*q/c;
    wtd_p=sorted_pvals./thresh;
    max_id=find(sorted_pvals<=thresh, 1, 'last' );
else
    error('Argument "method" needs to be "pdep" or "dep"');
end

if isempty(max_id)
    h=zeros(m,1);
    crit_p=0;
else
    h=pvals<=thresh(max_id);
    crit_p=sorted_pvals(max_id);
end

adj_p=zeros(m,1);
for i=m:-1:1
    adj_p(i)=min(1,min(wtd_p(i:end)));
end
adj_p=adj_p(unsort_order);

if strcmpi(report,'yes')
    n_sig=sum(h);
    if n_sig==1
        fprintf('Out of %d tests, %d is significant using the FDR at q = %.4f\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using the FDR at q = %.4f\n',m,n_sig,q);
    end
end
adj_ci_cvrg=[]; % Not implemented

end
