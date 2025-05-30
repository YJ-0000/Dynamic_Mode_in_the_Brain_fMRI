clear; clc;
%% load data
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_subExclude_indiv_10_B

[sub_ids,sorted_idx] = sort(sub_ids);
D = D(:,sorted_idx);
B_mean = B_mean(:,sorted_idx);
B_mean = B_mean ./ sum(B_mean,1);
B_var = B_var(:,sorted_idx);

D(1,:) = [];
D(:,tau==0) = [];
B_mean(:,tau==0) = [];
B_var(:,tau==0) = [];
sub_ids(tau==0) = [];

% D(:,tau<2000) = [];
% B_mean(:,tau<2000) = [];
% B_var(:,tau<2000) = [];
% sub_ids(tau<2000) = [];
%%
behav_factor = readtable('./data/scores_04_selected.csv');
% behav_factor = readtable('./data/scores_05_selected.csv');

%%
load secure_data/path_info;
gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
freesurfer_data_table = readtable(freesurfer_data_path,'VariableNamingRule','preserve');
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
for nrow = size(behav_factor,1):-1:1
    if ~any(sub_ids==behav_factor(nrow,'Subject').Variables)
        behav_factor(nrow,:) = [];
    end
end

for nrow = size(freesurfer_data_table,1):-1:1
    if ~any(sub_ids==freesurfer_data_table(nrow,'Subject').Variables)
        freesurfer_data_table(nrow,:) = [];
    end
end
gene_data_table = sortrows(gene_data_table, 'Subject');
behav_data_table = sortrows(behav_data_table, 'Subject');
behav_factor = sortrows(behav_factor, 'Subject');
freesurfer_data_table = sortrows(freesurfer_data_table, 'Subject');

age = gene_data_table.Age_in_Yrs;
sex = strcmp(behav_data_table.Gender,'F');
bmi = gene_data_table.BMI;
race1 = strcmp(gene_data_table.Race,'Am. Indian/Alaskan Nat.');
race2 = strcmp(gene_data_table.Race,'Asian/Nat. Hawaiian/Othr Pacific Is.');
race3 = strcmp(gene_data_table.Race,'Black or African Am.');
race4 = strcmp(gene_data_table.Race,'More than one');
race5 = strcmp(gene_data_table.Race,'Unknown or Not Reported');
race6 = strcmp(gene_data_table.Race,'White');
ICV = freesurfer_data_table.FS_IntraCranial_Vol;
TGMV = freesurfer_data_table.FS_Total_GM_Vol;
handedness = gene_data_table.Handedness;
edu_year = gene_data_table.SSAGA_Educ;



X_rest = [ones(size(age,1),1),age,sex, ICV.^(1/3),TGMV.^(1/3),handedness,tau(tau~=0)>2000];
% X_rest = [ones(size(age,1),1),age,sex, ICV.^(1/3),TGMV.^(1/3),handedness];
% X_rest = [ones(size(age,1),1),age,sex, ICV.^(1/3),TGMV.^(1/3)];
% X_rest = [ones(size(age,1),1),age,sex];
% X_rest = [ones(size(age,1),1)];


nan_sub = sum(isnan(X_rest),2)>0;
X_rest(nan_sub,:) = [];
gene_data_table(nan_sub,:) = [];
D(:,nan_sub) = [];
B_mean(:,nan_sub) = [];
B_var(:,nan_sub) = [];
behav_data_table(nan_sub,:) = [];
behav_factor(nan_sub,:) = [];

abs_D = abs(D(1:2:end,:))';
angle_D = angle(D(1:2:end,:))';
BB_mean = B_mean(1:2:end,:)';

idx_outlier = zeros(size(BB_mean,1),1);
% for n_dm = 1:5
%     idx_outlier = idx_outlier + isoutlier(BB_mean(:,n_dm));
%     idx_outlier = idx_outlier + isoutlier(abs_D(:,n_dm));
%     idx_outlier = idx_outlier + isoutlier(angle_D(:,n_dm));
% end

%% 
num_FA = size(behav_factor,2)-1;
num_DM = size(abs_D,2);
t_values_magnitude = zeros(num_DM, num_FA);
t_values_angle = zeros(num_DM, num_FA);
t_values_B_mean = zeros(num_DM, num_FA);
p_values_magnitude = zeros(num_DM, num_FA);
p_values_angle = zeros(num_DM, num_FA);
p_values_B_mean = zeros(num_DM, num_FA);

r_values_magnitude = zeros(num_DM, num_FA);
r_values_angle = zeros(num_DM, num_FA);
r_values_B_mean = zeros(num_DM, num_FA);

FA_values = table2array(behav_factor(:,2:end));

for n_f = 1:num_FA
    idx_outlier = idx_outlier + isoutlier(FA_values(:,n_f));
end
idx_outlier = idx_outlier > 0;

BB_mean(idx_outlier,:) = [];
abs_D(idx_outlier,:) = [];
angle_D(idx_outlier,:) = [];
FA_values(idx_outlier,:) = [];
X_rest(idx_outlier,:) = [];

BB_mean_resid = BB_mean - X_rest * pinv(X_rest) * BB_mean;
abs_D_resid = abs_D - X_rest * pinv(X_rest) * abs_D;
angle_D_resid = angle_D - X_rest * pinv(X_rest) * angle_D;

for n_dm = 1:num_DM
    for n_fa = 1:num_FA
        % Fit the GLM model
        current_predictors = [X_rest, BB_mean_resid(:, n_dm)];
        [b, d, st] = glmfit(current_predictors, FA_values(:, n_fa),'normal','Constant','off');
        t_B_mean = st.t(end);
        p_B_mean = st.p(end);
        FA_resid = FA_values(:, n_fa) - X_rest * pinv(X_rest) * FA_values(:, n_fa);
        r_B_mean = corrcoef(FA_resid, BB_mean(:, n_dm));
        
        current_predictors = [X_rest, abs_D_resid(:, n_dm)];
        [b, d, st] = glmfit(current_predictors, FA_values(:, n_fa),'normal','Constant','off');
        t_magnitude = st.t(end);
        p_magnitude = st.p(end);
        r_magnitude = corrcoef(FA_resid, abs_D(:, n_dm));

        current_predictors = [X_rest, angle_D_resid(:, n_dm)];
        [b, d, st] = glmfit(current_predictors, FA_values(:, n_fa),'normal','Constant','off');
        t_angle = st.t(end);
        p_angle = st.p(end);
        r_angle = corrcoef(FA_resid, angle_D(:, n_dm));
        
        t_values_magnitude(n_dm, n_fa) = t_magnitude;
        t_values_angle(n_dm, n_fa) = t_angle;
        t_values_B_mean(n_dm, n_fa) = t_B_mean;
        p_values_magnitude(n_dm, n_fa) = p_magnitude;
        p_values_angle(n_dm, n_fa) = p_angle;
        p_values_B_mean(n_dm, n_fa) = p_B_mean;
        r_values_magnitude(n_dm, n_fa) = r_magnitude(2);
        r_values_angle(n_dm, n_fa) = r_angle(2);
        r_values_B_mean(n_dm, n_fa) = r_B_mean(2);
    end
end
DM_names = {'Principal','SN-to-DMN','Bi-asym','FV-to-SM','SN-to-CEN'};
if num_FA == 4
    col_names = {'Mental Health', 'Cognition', 'Processing Speed', 'Substance Use'};
elseif num_FA == 5
    col_names = {'Mental Health', 'Cognition', 'Internalizing','Processing Speed', 'Substance Use'};
end

disp(' ');
fprintf('############ Mode-engagement level #############\n')
disp(array2table(t_values_B_mean, 'VariableNames', col_names, 'RowNames', DM_names));
fprintf('############ Persistence rate #############\n')
disp(array2table(t_values_magnitude, 'VariableNames', col_names, 'RowNames', DM_names));
fprintf('############ Progression rate #############\n')
disp(array2table(t_values_angle, 'VariableNames', col_names, 'RowNames', DM_names));

%%
p_concat = [p_values_B_mean(:); p_values_magnitude(:); p_values_angle(:)];
FDR = mafdr(p_concat);
FDR_B_mean = reshape(FDR(1:(num_FA*num_DM)),[num_DM,num_FA]);
FDR_magnitude = reshape(FDR((num_FA*num_DM+1):(2*num_FA*num_DM)),[num_DM,num_FA]);
FDR_angle = reshape(FDR((2*num_FA*num_DM+1):(3*num_FA*num_DM)),[num_DM,num_FA]);
fprintf('############ Mode-engagement level #############\n')
disp(array2table(FDR_B_mean, 'VariableNames', col_names, 'RowNames', DM_names));
fprintf('############ Persistence rate #############\n')
disp(array2table(FDR_magnitude, 'VariableNames', col_names, 'RowNames', DM_names));
fprintf('############ Progression rate #############\n')
disp(array2table(FDR_angle, 'VariableNames', col_names, 'RowNames', DM_names));
%%
[II, JJ] = ndgrid(1:num_DM, 1:num_FA);
II = II(:);
JJ = JJ(:);
num_iterations = length(II);

num_sub = size(BB_mean,1);

num_perm = 10000;
max_t_list = zeros(1,num_perm);
for n_perm = 1:num_perm
    fprintf('%d ', n_perm)
    t_values_magnitude_perm = zeros(num_DM, num_FA);
    t_values_angle_perm = zeros(num_DM, num_FA);
    t_values_B_mean_perm = zeros(num_DM, num_FA);
%     t_values_B_var_perm = zeros(num_abs, num_FA);
    for kk = 1:num_iterations
        ii = II(kk);
        jj = JJ(kk);
        
        perm_idx = randperm(num_sub);

        % Fit the GLM model
        current_predictors = [X_rest, abs_D_resid(perm_idx, ii)];
        [b, d, st] = glmfit(current_predictors, FA_values(:, jj),'normal','Constant','off');
        t_magnitude = st.t(end);

        current_predictors = [X_rest, angle_D_resid(perm_idx, ii)];
        [b, d, st] = glmfit(current_predictors, FA_values(:, jj),'normal','Constant','off');
        t_angle = st.t(end);
        
        current_predictors = [X_rest, BB_mean_resid(perm_idx, ii)];
        [b, d, st] = glmfit(current_predictors, FA_values(:, jj),'normal','Constant','off');
        t_B_mean = st.t(end);
        
%         current_predictors = [X_rest, BB_var(perm_idx, ii)];
%         [b, d, st] = glmfit(current_predictors, FA_values(:, jj),'normal','Constant','off');
%         t_B_var = st.t(end);

        % Store the p-values in preallocated matrices
        t_values_magnitude_perm(ii, jj) = t_magnitude;
        t_values_angle_perm(ii, jj) = t_angle;
        t_values_B_mean_perm(ii, jj) = t_B_mean;
%         t_values_B_var_perm(ii,jj) = t_B_var;
    end
    max_t_list(n_perm) = max(abs([t_values_magnitude_perm(:);t_values_angle_perm(:);t_values_B_mean_perm(:)])); 
end

p_values_magnitude_FWE = zeros(num_DM, num_FA);
p_values_angle_FWE = zeros(num_DM, num_FA);
p_values_B_mean_FWE = zeros(num_DM, num_FA);
% p_values_B_var_FWE = zeros(num_abs, num_FA);
for ii = 1:num_DM
    for jj = 1:num_FA
        p_values_magnitude_FWE(ii, jj) = mean(max_t_list>abs(t_values_magnitude(ii, jj)));
        p_values_angle_FWE(ii, jj) = mean(max_t_list>abs(t_values_angle(ii, jj)));
        p_values_B_mean_FWE(ii, jj) = mean(max_t_list>abs(t_values_B_mean(ii, jj)));
%         p_values_B_var_FWE(ii, jj) = mean(max_t_list>abs(t_values_B_var(ii, jj)));
    end
end
disp(' ');
fprintf('############ Mode-engagement level #############\n')
disp(array2table(p_values_B_mean_FWE, 'VariableNames', col_names, 'RowNames', DM_names));
fprintf('############ Persistence rate #############\n')
disp(array2table(p_values_magnitude_FWE, 'VariableNames', col_names, 'RowNames', DM_names));
fprintf('############ Progression rate #############\n')
disp(array2table(p_values_angle_FWE, 'VariableNames', col_names, 'RowNames', DM_names));


%% Display (heatmap)
DM_order = [1,2,5,4,3];
DM_names_correct_order = DM_names(DM_order);

temp_max = max(abs([t_values_B_mean(:); t_values_magnitude(:); t_values_angle(:)]));
minValue = -temp_max;
maxValue = temp_max;

for ii = 1:3

    % Create a new figure window
    figure;
    
    if ii == 1
        t_values_mat = t_values_B_mean;
        p_FWE_mat = p_values_B_mean_FWE;
        q_FDR_mat = FDR_B_mean;
    elseif ii ==2
        t_values_mat = t_values_magnitude;
        p_FWE_mat = p_values_magnitude_FWE;
        q_FDR_mat = FDR_magnitude;
    else
        t_values_mat = t_values_angle;
        p_FWE_mat = p_values_angle_FWE;
        q_FDR_mat = FDR_angle;
    end
    
    t_values_mat = t_values_mat(DM_order,:);
    p_FWE_mat = p_FWE_mat(DM_order,:);
    q_FDR_mat = q_FDR_mat(DM_order,:);

    % Display the matrix as an image
    imagesc(t_values_mat);
    caxis([minValue maxValue]);
    
    % Positions (normalized between 0 and 1)
    positions = [0, 0.5, 1];
    % Define the number of colors in the colormap
    nColors = 256;
    % Corresponding RGB colors
    % Blue: [0.2298, 0.2980, 0.7530]
    % White: [1, 1, 1]
    % Red: [0.7530, 0.2314, 0.0980]
    colors = [
        0.2298, 0.2980, 0.7530;    % Blue
        1.0000, 1.0000, 1.0000;    % White
        0.7530, 0.1504, 0.0980     % Red
    ];

    % Create a matrix for interpolation
    % Each row represents R, G, B channels respectively
    R = colors(:,1);
    G = colors(:,2);
    B = colors(:,3);

    % Define the query points for interpolation
    x = linspace(0, 1, nColors);

    % Interpolate each color channel
    R_interp = interp1(positions, R, x, 'linear');
    G_interp = interp1(positions, G, x, 'linear');
    B_interp = interp1(positions, B, x, 'linear');

    % Combine the interpolated channels into a colormap
    coolwarm = [R_interp', G_interp', B_interp'];

    % Apply the custom 'coolwarm' colormap
    colormap(coolwarm);
    
    save results/colormap_coolwarm coolwarm

    % Add a colorbar to the figure
    colorbar;

    % Get the current axes handle
    ax = gca;

    % Set the font name and size for the axes
    ax.FontName = 'Times New Roman';
    ax.FontSize = 24;

    % Define X-axis ticks and labels
    ax.XTick = 1:numel(col_names);
    ax.XTickLabel = col_names;

    % Optionally, rotate X-axis labels for better readability
    % Uncomment the next line if you want to rotate the labels by 45 degrees
    % ax.XTickLabelRotation = 45;

    % Define Y-axis ticks and labels
    ax.YTick = 1:numel(DM_names_correct_order);
    ax.YTickLabel = DM_names_correct_order;

    % Hold the current plot to add text annotations
    hold on;
    
    y_offset = 0.1;
    % Loop through each element in the t_values_B_mean matrix
    for i = 1:size(p_FWE_mat, 1)
        for j = 1:size(p_FWE_mat, 2)
            % Check if the corresponding p-value is less than 0.05
            if p_FWE_mat(i,j) < 0.001
                text(j, i + y_offset, '***', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 20, ...                    
                    'Color', 'k');                         % Set text color to red
            elseif p_FWE_mat(i,j) < 0.01
                text(j, i + y_offset, '**', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 20, ...                    
                    'Color', 'k');                         % Set text color to red
            elseif p_FWE_mat(i,j) < 0.05
                text(j, i + y_offset, '*', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 20, ...                    
                    'Color', 'k');                         % Set text color to red
            elseif q_FDR_mat(i,j) < 0.05
                text(j, i + y_offset, '+', ...
                    'HorizontalAlignment', 'center', ...
                    'FontSize', 20, ...                    
                    'Color', 'k');                         % Set text color to red
            end
        end
    end

    % Release the hold on the current plot
    hold off;

    % Optional: Adjust the figure's overall properties for better aesthetics
    set(gcf, 'Color', 'w'); % Set figure background to white
end

%% scatter plots

DM_order = [1,2,5,4,3];
for n_metric = 1:3
    if n_metric == 1
        metric = BB_mean;
        metric_resid = BB_mean_resid;
        metric_name = 'BB\_mean';
    elseif n_metric == 2
        metric = abs_D;
        metric_resid = abs_D_resid;
        metric_name = 'abs\_D';
    else
        metric = angle_D;
        metric_resid = angle_D_resid;
        metric_name = 'angle\_D';
    end
    
    figure('Name', metric_name);
    plot_n = 0;
    for n_dm = DM_order
        for n_fa = 1:num_FA
            plot_n = plot_n + 1;
            ax = subplot(num_DM, num_FA, plot_n);
            
            % 잔차 계산
            FA_resid = FA_values(:, n_fa) - X_rest * pinv(X_rest) * FA_values(:, n_fa);
%             x = metric(:, n_dm);
%             y = FA_resid;
            x = FA_resid;
            y = metric(:, n_dm);
            
            % 산점도
            scatter(x, y, 'filled');
            hold on;
            
            % 선형 회귀선 계산 및 그리기
            p = polyfit(x, y, 1);       % 1차 회귀
            yfit = polyval(p, x);
            plot(x, yfit, 'r-', 'LineWidth', 1.5);
            
            % 상관계수 계산
            R = corrcoef(x, y);
            r = R(1,2);
            
            % 축 범위 가져와서 텍스트 위치 계산
            xl = xlim(ax);
            yl = ylim(ax);
            xpos = xl(1) + 0.05 * diff(xl);
            ypos = yl(2) - 0.10 * diff(yl);
            text(xpos, ypos, sprintf('r = %.2f', r), 'FontSize', 10, 'Color', 'r');
            
            % 레이블
%             xlabel(sprintf('Metric %d (DM=%d)', n_metric, n_dm));
%             ylabel(sprintf('FA resid (FA=%d)', n_fa));
%             title(sprintf('FA %d vs %s (DM %d)', n_fa, metric_name, n_dm));
            
            hold off;
        end
    end
end

