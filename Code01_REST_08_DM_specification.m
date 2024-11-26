clear; clc; close all;
%%
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10_B
%%
load results/RSFC_standard_grad
load results/RSFC

%% FC
N = 360;

max_dm = 5;

Phi_DM_ang_accum = zeros(N,N,max_dm);

for n_dm = 1:max_dm %1:5
        
    Phi_DM_ang = angle(Phi_sorted(1:N,2*n_dm));
    
    Phi_DM_ang_mat = cos(Phi_DM_ang - Phi_DM_ang');
    
    Phi_DM_ang_accum(:,:,n_dm) = Phi_DM_ang_mat;
    
    
    index_tril = tril(ones(size(Phi_DM_ang_mat)),-1);
    Phi_DM_ang_vector = Phi_DM_ang_mat(index_tril ~= 0);
    group_corr_mat_vector = group_corr_mat_original(index_tril ~= 0);
    
    Phi_DM_ang_accum_vecter = squeeze(mean(Phi_DM_ang_accum(:,:,1:n_dm),3));
    Phi_DM_ang_accum_vecter = Phi_DM_ang_accum_vecter(index_tril ~= 0);
    
    [r,p] = corrcoef(Phi_DM_ang_vector,group_corr_mat_vector);
    
    disp(['DM',num2str(n_dm),' - FC correlation with cos(ang diff): ',num2str(r(2),'%.6f')]);
    
    % Compute modularity and sort matrices
    [ci, Q] = modularity_und(group_corr_mat_original);
    [ci_recon, Q_recon] = modularity_und(Phi_DM_ang_mat);
%     [ci] = community_louvain(group_corr_mat_original);
%     [ci_recon] = community_louvain(Phi_DM_ang_mat);
    
    % Sort the reconstructed FC matrix
    [~, idx_sorted] = sort(ci);
    group_corr_mat_original_sorted = group_corr_mat_original(idx_sorted, idx_sorted);
    sorted_ci = ci(idx_sorted); % Sorted community indices

    % Sort the original group correlation matrix
    Phi_DM_ang_mats_sorted = Phi_DM_ang_mat(idx_sorted, idx_sorted);

    % Display community information
    disp(['The number of matching community identities: ', num2str(sum(ci == ci_recon))]);
    disp(['The NMI of matching community identities: ', num2str(nmi(ci, ci_recon))]);
    
    
    % Define red-white-blue colormap
    n = 256; % Number of colors
    half = floor(n/2);

    % Blue to white
    cmap1 = [linspace(0,1,half)', linspace(0,1,half)', linspace(1,1,half)'];

    % White to red
    cmap2 = [linspace(1,1,n-half)', linspace(1,0,n-half)', linspace(1,0,n-half)'];

    % Combine to create red-white-blue colormap
    redWhiteBlue = [cmap1; cmap2];

    % Create figure with specified size to ensure square heatmaps
    figure('Position', [100, 100, 1400, 700]); % [left, bottom, width, height]

    % Define index for display labels (assuming 360 nodes as per original code)
    total_nodes = size(group_corr_mat_original, 1); % Generalizing to any number of nodes
    idx_labels = rem(1:total_nodes,100) == 0;

    %%% Plot the first heatmap (Reconstructed FC) using imagesc
    subplot(1,2,1); 
    imagesc(group_corr_mat_original_sorted);
    colormap(redWhiteBlue);
    colorbar;
    axis square;
    hold on; % Enable adding lines

    % Set labels with desired font
    xlabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    ylabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    % title('Reconstructed FC', 'FontName', 'Times New Roman', 'FontSize', 24);

    % Customize tick labels
    set(gca, 'XTick', 1:total_nodes, 'YTick', 1:total_nodes, ...
             'FontName', 'Times New Roman', 'FontSize', 24);
    set(gca, 'XTickLabel', repmat({''}, 1, total_nodes));
    set(gca, 'YTickLabel', repmat({''}, 1, total_nodes));
    if any(idx_labels)
        xticks_labels = find(idx_labels);
        yticks_labels = find(idx_labels);
        set(gca, 'XTick', xticks_labels, 'YTick', yticks_labels);
        set(gca, 'XTickLabel', num2str(xticks_labels'));
        set(gca, 'YTickLabel', num2str(yticks_labels'));
    end

    %%% Plot the second heatmap (Original Group Correlation) using imagesc
    subplot(1,2,2); 
    imagesc(Phi_DM_ang_mats_sorted);
    colormap(redWhiteBlue);
    colorbar;
    axis square;
    hold on; % Enable adding lines

    % Set labels with desired font
    xlabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    ylabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    % title('Original Group Correlation', 'FontName', 'Times New Roman', 'FontSize', 24);

    % Customize tick labels
    set(gca, 'XTick', 1:total_nodes, 'YTick', 1:total_nodes, ...
             'FontName', 'Times New Roman', 'FontSize', 24);
    set(gca, 'XTickLabel', repmat({''}, 1, total_nodes));
    set(gca, 'YTickLabel', repmat({''}, 1, total_nodes));
    if any(idx_labels)
        xticks_labels = find(idx_labels);
        yticks_labels = find(idx_labels);
        set(gca, 'XTick', xticks_labels, 'YTick', yticks_labels);
        set(gca, 'XTickLabel', num2str(xticks_labels'));
        set(gca, 'YTickLabel', num2str(yticks_labels'));
    end

    % Optional: Adjust layout for better appearance
    sgtitle('Comparison of Original FC and angle difference map', ...
        'FontName', 'Times New Roman', 'FontSize', 24);
    
end


%% FC accumulation
lambda_indiv_mean = mean(D(2:end,:),2);
[lambda_indiv_mean_sorted,idx_sorted] = sort(lambda_indiv_mean,'descend');
Phi_sorted_indiv = Phi_sorted(:,idx_sorted);

B_group = mean(B_mean(:,tau~=0),2);

max_dm = 5;

frame_dt = 0.5;
frame_num = 10000;
TRtarget = 1.5;
N = 360;
reconstructed_timeseries = zeros(N,frame_num+1,max_dm);
for n_dm = 1:max_dm
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num)/abs(lambda(DM_conjugate1_num));
    lambda_conjugate2 = lambda(DM_conjugate2_num)/abs(lambda(DM_conjugate2_num));
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
%         reconstructed_timeseries(:,frame+1,n_dm) = abs(lambda(DM_conjugate1_num))*eigenstate;
        reconstructed_timeseries(:,frame+1,n_dm) = B_group(DM_conjugate1_num)*eigenstate;
    end
end

r_accum = zeros(1,max_dm);
for n_dm = 1:max_dm
    reconstructed_FC = corrcoef(squeeze(mean(reconstructed_timeseries(:,:,1:n_dm),3))');
    index_tril = tril(ones(size(reconstructed_FC)),-1);
    reconstructed_FC_vector = reconstructed_FC(index_tril ~= 0);
    group_corr_mat_vector = group_corr_mat_original(index_tril ~= 0);
%     group_corr_mat_vector = group_corr_mat_global_removed(index_tril ~= 0);

    r = corrcoef(reconstructed_FC_vector,group_corr_mat_vector);
    r_accum(n_dm) = r(2);
end

%% calculate dFC (DM1)
N_cortex = 360;
WL = round(30*0.72/0.5);
dFC_reconstructed_DM1 = zeros(N_cortex,N_cortex);
for roi1 = 1:N_cortex
    fprintf('%d ',roi1);
    x = squeeze(reconstructed_timeseries(roi1,:,1));
    parfor roi2 = (roi1+1):N_cortex
        y = squeeze(reconstructed_timeseries(roi2,:,1));
        dFC_reconstructed_DM1(roi1,roi2) = compute_dstd(x, y, WL);
    end
end
fprintf('\n');

dFC_reconstructed_DM1 = dFC_reconstructed_DM1 + dFC_reconstructed_DM1';

save DMs/temp_dFC_DMs dFC_reconstructed_DM1

%% calculate dFC
N_cortex = 360;
WL = round(30*0.72/0.5);
dFC_reconstructed = zeros(N_cortex,N_cortex);
for roi1 = 1:N_cortex
    fprintf('%d ',roi1);
    x = squeeze(mean(reconstructed_timeseries(roi1,:,1:max_dm),3));
    parfor roi2 = (roi1+1):N_cortex
        y = squeeze(mean(reconstructed_timeseries(roi2,:,1:max_dm),3));
        dFC_reconstructed(roi1,roi2) = compute_dstd(x, y, WL);
    end
end
fprintf('\n');

dFC_reconstructed = dFC_reconstructed + dFC_reconstructed';

save('DMs/temp_dFC_DMs', 'dFC_reconstructed','-append');

%% Community analysis
load results/RSFC

disp(['The correlation between reconstructed FC and FC = ', num2str(r(2),'%0.4f')]);
figure; plot(r_accum); title('recon FC vs FC'); xlabel('number of DMs');

% figure; scatter(reconstructed_FC_vector,group_corr_mat_vector);

% Compute modularity and sort matrices
[ci, Q] = modularity_und(group_corr_mat_original);
[ci_recon, Q_recon] = modularity_und(reconstructed_FC);
ci_recon2 = ci_recon;
ci_recon(ci_recon2==1) = 2;
ci_recon(ci_recon2==2) = 1;
ci_recon(ci_recon2==3) = 3;

N = 360;
label_idx_list = 1:N;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
ci_grad = zeros(size(labels.cdata));
for n_reg = 1:N
    ci_grad(labels.cdata==label_idx_list(n_reg)) = ci(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, ci_grad, 'dscalar');
cifti_write(lag1, ['results/community_structure.dscalar.nii']);

% Sort the reconstructed FC matrix
[~, idx_sorted] = sort(ci);
group_corr_mat_original_sorted = group_corr_mat_original(idx_sorted, idx_sorted);
sorted_ci = ci(idx_sorted); % Sorted community indices

% Sort the original group correlation matrix
[~, idx_sorted_recon] = sort(ci_recon);
reconstructed_FC_sorted = reconstructed_FC(idx_sorted_recon, idx_sorted_recon);
sorted_ci_recon = ci_recon(idx_sorted_recon); % Sorted community indices for original

% Display community information
disp(['The number of matching community identities: ', num2str(sum(ci == ci_recon))]);
disp(['The NMI of matching community identities: ', num2str(nmi(ci, ci_recon))]);

% Define red-white-blue colormap
n = 256; % Number of colors
half = floor(n/2);

% Blue to white
cmap1 = [linspace(0,1,half)', linspace(0,1,half)', linspace(1,1,half)'];

% White to red
cmap2 = [linspace(1,1,n-half)', linspace(1,0,n-half)', linspace(1,0,n-half)'];

% Combine to create red-white-blue colormap
redWhiteBlue = [cmap1; cmap2];

% Create figure with specified size to ensure square heatmaps
figure('Position', [100, 100, 1400, 700]); % [left, bottom, width, height]

% Define index for display labels (assuming 360 nodes as per original code)
total_nodes = size(reconstructed_FC, 1); % Generalizing to any number of nodes
idx_labels = rem(1:total_nodes,100) == 0;

%%% Plot the first heatmap (Reconstructed FC) using imagesc
subplot(1,2,1); 
imagesc(group_corr_mat_original_sorted);
colormap(redWhiteBlue);
colorbar;
axis square;
hold on; % Enable adding lines

% Set labels with desired font
xlabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
% title('Reconstructed FC', 'FontName', 'Times New Roman', 'FontSize', 24);

% Customize tick labels
set(gca, 'XTick', 1:total_nodes, 'YTick', 1:total_nodes, ...
         'FontName', 'Times New Roman', 'FontSize', 24);
set(gca, 'XTickLabel', repmat({''}, 1, total_nodes));
set(gca, 'YTickLabel', repmat({''}, 1, total_nodes));
if any(idx_labels)
    xticks_labels = find(idx_labels);
    yticks_labels = find(idx_labels);
    set(gca, 'XTick', xticks_labels, 'YTick', yticks_labels);
    set(gca, 'XTickLabel', num2str(xticks_labels'));
    set(gca, 'YTickLabel', num2str(yticks_labels'));
end

% Calculate cluster boundaries for the first heatmap
clusters = unique(sorted_ci);
cluster_sizes = arrayfun(@(c) sum(sorted_ci == c), clusters);
boundary_positions = cumsum(cluster_sizes);

% Draw boundary lines
start_pos = 1;
for i = 1:length(clusters)
    cluster_size = cluster_sizes(i);
    % Define rectangle position [x, y, width, height]
    rect = [start_pos - 0.5, start_pos - 0.5, cluster_size, cluster_size];
    % Draw rectangle
    rectangle('Position', rect, 'EdgeColor', 'k', 'LineWidth', 1);
    % Update starting position
    start_pos = start_pos + cluster_size;
end

hold off; % Disable adding more lines

%%% Plot the second heatmap (Original Group Correlation) using imagesc
subplot(1,2,2); 
imagesc(reconstructed_FC_sorted);
colormap(redWhiteBlue);
colorbar;
axis square;
hold on; % Enable adding lines

% Set labels with desired font
xlabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
% title('Original Group Correlation', 'FontName', 'Times New Roman', 'FontSize', 24);

% Customize tick labels
set(gca, 'XTick', 1:total_nodes, 'YTick', 1:total_nodes, ...
         'FontName', 'Times New Roman', 'FontSize', 24);
set(gca, 'XTickLabel', repmat({''}, 1, total_nodes));
set(gca, 'YTickLabel', repmat({''}, 1, total_nodes));
if any(idx_labels)
    xticks_labels = find(idx_labels);
    yticks_labels = find(idx_labels);
    set(gca, 'XTick', xticks_labels, 'YTick', yticks_labels);
    set(gca, 'XTickLabel', num2str(xticks_labels'));
    set(gca, 'YTickLabel', num2str(yticks_labels'));
end

% Calculate cluster boundaries for the second heatmap
clusters_recon = unique(sorted_ci_recon);
cluster_sizes_recon = arrayfun(@(c) sum(sorted_ci_recon == c), clusters_recon);
boundary_positions_recon = cumsum(cluster_sizes_recon);

% Draw boundary lines
start_pos = 1;
for i = 1:length(clusters)
    cluster_size = cluster_sizes_recon(i);
    % Define rectangle position [x, y, width, height]
    rect = [start_pos - 0.5, start_pos - 0.5, cluster_size, cluster_size];
    % Draw rectangle
    rectangle('Position', rect, 'EdgeColor', 'k', 'LineWidth', 1);
    % Update starting position
    start_pos = start_pos + cluster_size;
end

hold off; % Disable adding more lines

% Optional: Adjust layout for better appearance
sgtitle('Comparison of Reconstructed FC and Original Group Correlation', ...
    'FontName', 'Times New Roman', 'FontSize', 24);




%% FC gradient
N = 360;

max_dm = 5;

Phi_DM_ang_accum = zeros(N,N,max_dm);
r_accum = zeros(1,max_dm);

for n_dm = 1:max_dm
        
    Phi_DM_ang = angle(Phi_sorted(1:N,2*n_dm));
    
    Phi_DM_ang_mat = cos(Phi_DM_ang - Phi_DM_ang');
    coeff = pca(Phi_DM_ang_mat);
    
    r = corrcoef(coeff(:,1),-RSFC_grad(:,1));
    
    disp(['Correlation between the real and reconstructed FC gradient of DM #',num2str(n_dm),': ',num2str(r(2),'%.4f')]);
end
%% lag projection
disp('###### Lag projection analysis #######');
 
N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

load('results/time_lags.mat');

label_idx_list = 1:N;
lag_data = mean(grp_lags_mean,1)';
grad1_data = zeros(size(labels.cdata));
for n_reg = 1:N
    grad1_data(labels.cdata==label_idx_list(n_reg)) = lag_data(n_reg);
end
lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
cifti_write(lag1, ['results/original_lag.dscalar.nii']);

for n_dm = 1:1 %1:5
    Phi_DM_mag = abs(Phi_sorted(1:N,2*n_dm));
        
    Phi_DM_ang = angle(Phi_sorted(1:N,2*n_dm));
    
    Phi_DM_ang_mat = Phi_DM_ang - Phi_DM_ang';
    
    lag_DM = mean(Phi_DM_ang_mat,2);
    
    lag_data = mean(grp_lags_mean,1)';
    
    [r,p] = corrcoef(lag_DM,lag_data);
    
    disp(['DM',num2str(n_dm),' - correlation of lag structure: ',num2str(r(2),'%.6f')]);
    
%     figure; scatter(lag_DM,lag_data);
    
    grad1_data = zeros(size(labels.cdata));
    for n_reg = 1:N
        grad1_data(labels.cdata==label_idx_list(n_reg)) = lag_DM(n_reg);
    end
    lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
    cifti_write(lag1, ['results/DM',num2str(n_dm),'_lag.dscalar.nii']);
    
    lag_DM_list(:,n_dm) = lag_DM;
end

%% global QPP
load results/QPP_100_noGSR
% lambda_indiv_mean = mean(D(2:end,:),2);

ref_t = 28;
frame_dt = 0.72;
frame_num = 278;
TRtarget = 1.5;
N = 360;
reconstructed_timeseries = zeros(N,frame_num+1);
for n_dm = 1:1
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);

    
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        reconstructed_timeseries(:,frame+1) = eigenstate;
    end
    
    corr_time_series_PrincipalGradient_DMN = zeros(1,frame_num+1);
    corr_time_series = zeros(1,frame_num-size(final_qpp_template,2));
    for frame = 1:size(corr_time_series,2)
        temp_dm = reconstructed_timeseries(:,frame:frame+size(final_qpp_template,2)-1);
        r = corrcoef(temp_dm(:),final_qpp_template(:));
        corr_time_series(frame) = r(2);
        r = corrcoef(reconstructed_timeseries(:,frame),RSFC_grad(:,1));
        corr_time_series_PrincipalGradient_DMN(frame+1) = -r(2);
    end
    [max_r,max_idx] = max(corr_time_series);
    disp(['Max correlation with QPP - DM',num2str(n_dm),': ',num2str(max_r,'%.4f'), ' at time t=', num2str(max_idx*frame_dt)]);
    figure; plot(frame_dt*(0:1:frame_num),corr_time_series_PrincipalGradient_DMN,'LineWidth',2);
    xline(max_idx*frame_dt,'r--');
    xline(max_idx*frame_dt+30*frame_dt,'r--');
end

%% anti-corr QPP
load results/QPP_100_GSR
lambda_indiv_mean = mean(D(2:end,:),2);
% B_group = mean(B_mean(:,tau~=0),2);

frame_dt = 0.72;
frame_num = 3000;
TRtarget = 1.5;
N = 360;
reconstructed_timeseries = zeros(N,frame_num+1,9);
tracking_num = 0;
for n_dm =[3, 1]
    ref_t = 0;
    tracking_num = tracking_num + 1;
    
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num)/abs(lambda(DM_conjugate1_num));
    lambda_conjugate2 = lambda(DM_conjugate2_num)/abs(lambda(DM_conjugate2_num));
    
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        reconstructed_timeseries(:,frame+1,tracking_num) = eigenstate;
    end
    
    accum_timeseries = squeeze(sum(reconstructed_timeseries(:,:,1:tracking_num),3));
    
    corr_time_series = zeros(1,frame_num-size(final_qpp_template,2));
    for frame = 1:size(corr_time_series,2)
        temp_dm = accum_timeseries(:,frame:frame+size(final_qpp_template,2)-1);
        r = corrcoef(temp_dm(:),final_qpp_template(:));
        corr_time_series(frame) = r(2);
    end
    [max_r,max_idx] = max(corr_time_series);
    disp(['Max correlation with QPP - DM',num2str(n_dm),': ',num2str(max_r,'%.4f'), ' at time t=', num2str(max_idx*frame_dt)]);
end

%% RSFC grad correlation time course
frame_dt = 0.5;
frame_num = 200;
TRtarget = 1.5;
N = 360;

% ref_t = 0;
load results/RSFC_standard_grad
load results/RSFC

load('results/FPN_vector.mat');
load('results/SMLV_vector.mat');
load('results/Task_negative_vector.mat');
load('results/Right_aIns_seed_vector.mat');
load('results/DAN_vector.mat');
load('results/VAN_vector.mat');
load('results/LAN_vector.mat');
load('results/Lateralized_index.mat');
load('results/Group_CAPs.mat');
load('results/Yeo_7networks.mat');

left_hemi_idx = 181:360;
right_hemi_idx = 1:180;

for n_dm = 1:5
    if n_dm == 1
        %%% principal DM
        ref_t = 27;
    elseif n_dm == 3
        ref_t = 31;
    elseif n_dm == 5
        ref_t = 44.5;
    elseif n_dm == 2
        ref_t = 56.5;
    elseif n_dm == 4
        ref_t = 14;
    else
        ref_t = 0;
    end
    
    
    corr_time_series_PrincipalGradient_DMN = zeros(1,frame_num+1);
    corr_time_series_SecondGradient = zeros(1,frame_num+1);
    corr_time_series_DMN = zeros(1,frame_num+1);
    corr_time_series_D_attention = zeros(1,frame_num+1);
    corr_time_series_V_attention = zeros(1,frame_num+1);
    corr_time_series_executive = zeros(1,frame_num+1);
    corr_time_series_salience = zeros(1,frame_num+1);
    corr_time_series_SMLV = zeros(1,frame_num+1);
    corr_time_series_LAN = zeros(1,frame_num+1);
    corr_time_series_MLI = zeros(1,frame_num+1);
    
    corr_time_series_PrincipalGradient_left = zeros(1,frame_num+1);
    corr_time_series_PrincipalGradient_right = zeros(1,frame_num+1);
    corr_time_series_executive_Yeo7_left = zeros(1,frame_num+1);
    corr_time_series_executive_Yeo7_right = zeros(1,frame_num+1);
    
    corr_time_series_CAP = zeros(size(CAPs,1),frame_num+1);
    corr_time_series_Yeo_7network = zeros(7,frame_num+1);
    
    Phi_DM_mag = abs(Phi_sorted(1:N,2*n_dm));
    
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    phase_lambda = angle(lambda_conjugate1);
    period = 1.5 * 2 * pi/phase_lambda;
    fprintf('The period of DM #%d: (period) %.2fs, (freq) %.6fHz \n',n_dm,period,1/period);
    
    RSFC_grad_temp = RSFC_grad(1:N,:);
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        
        r = corrcoef(eigenstate,RSFC_grad_temp(:,1));
        corr_time_series_PrincipalGradient_DMN(frame+1) = -r(2);
        
        r = corrcoef(eigenstate,RSFC_grad_temp(:,2));
        corr_time_series_SecondGradient(frame+1) = r(2);
        
        r = corrcoef(eigenstate,Task_negative_vector);
        corr_time_series_DMN(frame+1) = r(2);
        
        r = corrcoef(eigenstate,DAN_vector);
        corr_time_series_D_attention(frame+1) = r(2);
        
        r = corrcoef(eigenstate,VAN_vector);
        corr_time_series_V_attention(frame+1) = r(2);
        
        r = corrcoef(eigenstate,FPN_vector);
        corr_time_series_executive(frame+1) = r(2);
        
        r = corrcoef(eigenstate,Salience_vector);
        corr_time_series_salience(frame+1) = r(2);
        
        r = corrcoef(eigenstate,SMLV_vector);
        corr_time_series_SMLV(frame+1) = r(2);
        
        r = corrcoef(eigenstate,LAN_vector);
        corr_time_series_LAN(frame+1) = r(2);
        
        r = corrcoef(eigenstate,MLI);
        corr_time_series_MLI(frame+1) = r(2);
        
        for n_cap = 1:size(CAPs,1)
            r = corrcoef(eigenstate,CAPs(n_cap,:)');
            corr_time_series_CAP(n_cap,frame+1) = r(2);
        end
        
        for n_cap = 1:7
            r = corrcoef(eigenstate,Yeo_7network_all_vector(:,n_cap));
            corr_time_series_Yeo_7network(n_cap,frame+1) = r(2);
        end
        
        r = corrcoef(eigenstate(left_hemi_idx),RSFC_grad_temp(left_hemi_idx,1));
        corr_time_series_PrincipalGradient_left(frame+1) = -r(2);
        r = corrcoef(eigenstate(right_hemi_idx),RSFC_grad_temp(right_hemi_idx,1));
        corr_time_series_PrincipalGradient_right(frame+1) = -r(2);
        r = corrcoef(eigenstate(left_hemi_idx),Yeo_7network_all_vector(left_hemi_idx,3));
        corr_time_series_executive_Yeo7_left(frame+1) = r(2);
        r = corrcoef(eigenstate(right_hemi_idx),Yeo_7network_all_vector(right_hemi_idx,3));
        corr_time_series_executive_Yeo7_right(frame+1) = r(2);
        
    end
    
    figure; hold on;
    plot(frame_dt*(0:1:frame_num),corr_time_series_PrincipalGradient_DMN,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_SecondGradient,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_DMN,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_executive,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_salience,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_LAN,'LineWidth',2);
%     plot(frame_dt*(0:1:frame_num),corr_time_series_SMLV,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_D_attention,'b--','LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_V_attention,'r--','LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_MLI,'g--','LineWidth',2);
    
    for n_cap = 1:size(CAPs,1)
        plot(frame_dt*(0:1:frame_num),corr_time_series_CAP(n_cap,:),'-.','LineWidth',1);
    end
    
    for n_cap = 1:7
        plot(frame_dt*(0:1:frame_num),corr_time_series_Yeo_7network(n_cap,:),':','LineWidth',2);
    end
    
%     plot(frame_dt*(0:1:frame_num),corr_time_series_PrincipalGradient_left,'.','LineWidth',2);
%     plot(frame_dt*(0:1:frame_num),corr_time_series_PrincipalGradient_right,'.','LineWidth',2);
%     plot(frame_dt*(0:1:frame_num),corr_time_series_executive_Yeo7_left,'.','LineWidth',2);
%     plot(frame_dt*(0:1:frame_num),corr_time_series_executive_Yeo7_right,'.','LineWidth',2);
    
%     legend('Principal Grad','Second Grad','PCC seed','FPN','salience','SMLV','D attention','V attention','Lateralized index');
    yline(0.7,'k--');


    legend('Principal Grad','Second Grad','DMN','FPN','salience','Language','D attention','V attention','Lateralized index',...
        'CAP1','CAP2','CAP3','CAP4','CAP5','CAP6','CAP7','CAP8',...
        'Yeo1-DMN1','Yeo2-SM','Yeo3-CEN','Yeo4-Vis','Yeo5-DAN','Yeo6-VAN','Yeo7-DMN2',...
        'Location','southeast');
    %         'Principal Grad (left)', 'Principal Grad (right)', 'Yeo3-CEN (left)','Yeo3-CEN (right)',...

    xlabel('time (t)');
    ylabel('correlation');
end

%% RSFC grad correlation time course (subcortical gradient)
frame_dt = 0.5;
frame_num = 200;
TRtarget = 1.5;
N = 716;

load results/RSFC_subcortical_grad

for n_dm = 1
    
    if n_dm == 1
        %%% principal DM
        ref_t = 27;
    else
        ref_t = 0;
    end
    
    corr_time_series_hippocampus = zeros(1,frame_num+1);
    corr_time_series_thalamus = zeros(1,frame_num+1);
    corr_time_series_striatum = zeros(1,frame_num+1);
    corr_time_series_cerebellum = zeros(1,frame_num+1);
    corr_time_series_amygdala = zeros(1,frame_num+1);
    corr_time_series_brainstem = zeros(1,frame_num+1);
    
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        
        n_grad = 1;
        
        r = corrcoef(eigenstate(target_roi_idx_hippocampus),RSFC_grad_hippocampus(:,n_grad));
        corr_time_series_hippocampus(frame+1) = r(2);
        r = corrcoef(eigenstate(target_roi_idx_thalamus),RSFC_grad_thalamus(:,n_grad));
        corr_time_series_thalamus(frame+1) = r(2);
        r = corrcoef(eigenstate(target_roi_idx_striatum),RSFC_grad_striatum(:,n_grad));
        corr_time_series_striatum(frame+1) = r(2);
        r = corrcoef(eigenstate(target_roi_idx_cerebellum),RSFC_grad_cerebellum(:,n_grad));
        corr_time_series_cerebellum(frame+1) = r(2);
        r = corrcoef(eigenstate(target_roi_idx_amygdala),RSFC_grad_amygdala(:,n_grad));
        corr_time_series_amygdala(frame+1) = r(2);
        r = corrcoef(eigenstate(target_roi_idx_brainstem),RSFC_grad_brainstem(:,n_grad));
        corr_time_series_brainstem(frame+1) = r(2);
        
    end
    
    figure; hold on;
    plot(frame_dt*(0:1:frame_num),corr_time_series_hippocampus,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_amygdala,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_thalamus,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_striatum,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_brainstem,'LineWidth',2);
    plot(frame_dt*(0:1:frame_num),corr_time_series_cerebellum,'LineWidth',2);
    
    result_agg = [corr_time_series_hippocampus',...
                  corr_time_series_amygdala',...
                  corr_time_series_thalamus',...
                  corr_time_series_striatum',...
                  corr_time_series_brainstem',...
                  corr_time_series_cerebellum'];
    
    yline(0.7,'k--');
    hold off;
    ylim([-1,1]);
    legend('Hippocampus','Thalamus','Striatum','Cerebellum','Amygdala','Brainstem');
end

%% dFC
load results/dstd_all.mat
load results/execursion_all.mat

frame_dt = 0.5;
frame_num = 200;
TRtarget = 1.5;
N = 360;
reconstructed_timeseries = zeros(N,frame_num+1);
for n_dm = 1:1
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
   
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);

    
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        reconstructed_timeseries(:,frame+1) = eigenstate;
    end
    
    dFC = zeros(N,N,frame_num-30);
    for frame = 1:size(dFC,3)
        r = corrcoef(reconstructed_timeseries(:,frame:frame+30)');
        dFC(:,:,frame) = r;
    end
    dFC_std = std(dFC,0,3);
    
    Phi_DM_ang = angle(Phi_sorted(1:360,2*n_dm));
    Phi_DM_ang_mat = abs(Phi_DM_ang - Phi_DM_ang');
    
    
    index_tril = tril(ones(size(Phi_DM_ang_mat)),-1);
    Phi_DM_ang_vector = Phi_DM_ang_mat(index_tril ~= 0);

    dFC = dstd_group_45;
%     dFC = execursion_group;
    dFC(logical(eye(size(dFC)))) = nan;

    dFC_vector = dFC(index_tril ~= 0);
    
    k = 3.587; x0 = pi/4; y0 = 0;
%     pred_dFC = 1 ./ (1 + exp(-k*(Phi_DM_ang_vector - x0))) + y0;
%     pred_dFC(Phi_DM_ang_vector>pi/2) = 1 ./ (1 + exp(-k*(pi-Phi_DM_ang_vector(Phi_DM_ang_vector>pi/2) - x0))) + y0;
    
    %%% for std
    pred_dFC = 0.52 * sin((Phi_DM_ang_mat-pi/4)/0.5)/4 + (0.52/2);
    %%% for excursion
%     pred_dFC = cos((Phi_DM_ang_mat));
    
    
    pred_dFC(logical(eye(size(dFC)))) = nan;
    pred_dFC_vector = pred_dFC(index_tril ~= 0);
    dFC_reconstructed_vector = dFC_reconstructed(index_tril ~= 0);
    dFC_reconstructed_DM1_vector = dFC_reconstructed_DM1(index_tril ~= 0);
    
    r = corrcoef(dFC_reconstructed_vector,dFC_vector);
    
    disp(['dFC calculated. The correlation: ',num2str(r(2))]);
    
    [ci, Q] = modularity_und(group_corr_mat_original);
    % Sort the reconstructed FC matrix
    [~, idx_sorted] = sort(ci);
    dFC = dFC(idx_sorted, idx_sorted);
    pred_dFC = pred_dFC(idx_sorted, idx_sorted);
    dFC_reconstructed(logical(eye(size(dFC_reconstructed)))) = nan;
    dFC_reconstructed_sorted = dFC_reconstructed(idx_sorted, idx_sorted);
    dFC_reconstructed_DM1_sorted = dFC_reconstructed_DM1(idx_sorted, idx_sorted);
    sorted_ci = ci(idx_sorted); % Sorted community indices

    
    % Sort the original group correlation matrix
    Phi_DM_ang_mats_sorted = Phi_DM_ang_mat(idx_sorted, idx_sorted);    
    
    % Define red-white-blue colormap
    n = 256; % Number of colors
    half = floor(n/2);

    grayMap = gray(256); % 256 levels of gray

    % Create figure with specified size to ensure square heatmaps
    figure('Position', [100, 100, 1400, 700]); % [left, bottom, width, height]

    % Define index for display labels (assuming 360 nodes as per original code)
    total_nodes = size(dFC, 1); % Generalizing to any number of nodes
    idx_labels = rem(1:total_nodes,100) == 0;

    %%% Plot the first heatmap (Reconstructed FC) using imagesc
    subplot(1,2,1); 
    imagesc(dFC);
    colormap(grayMap);
    colorbar;
    axis square;
    hold on; % Enable adding lines

    % Set labels with desired font
    xlabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    ylabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    % title('Reconstructed FC', 'FontName', 'Times New Roman', 'FontSize', 24);

    % Customize tick labels
    set(gca, 'XTick', 1:total_nodes, 'YTick', 1:total_nodes, ...
             'FontName', 'Times New Roman', 'FontSize', 24);
    set(gca, 'XTickLabel', repmat({''}, 1, total_nodes));
    set(gca, 'YTickLabel', repmat({''}, 1, total_nodes));
    if any(idx_labels)
        xticks_labels = find(idx_labels);
        yticks_labels = find(idx_labels);
        set(gca, 'XTick', xticks_labels, 'YTick', yticks_labels);
        set(gca, 'XTickLabel', num2str(xticks_labels'));
        set(gca, 'YTickLabel', num2str(yticks_labels'));
    end

    %%% Plot the second heatmap (Original Group Correlation) using imagesc
    subplot(1,2,2); 
    imagesc(dFC_reconstructed_sorted);
    colormap(grayMap);
    colorbar;
    axis square;
    hold on; % Enable adding lines

    % Set labels with desired font
    xlabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    ylabel('Nodes', 'FontName', 'Times New Roman', 'FontSize', 24);
    % title('Original Group Correlation', 'FontName', 'Times New Roman', 'FontSize', 24);

    % Customize tick labels
    set(gca, 'XTick', 1:total_nodes, 'YTick', 1:total_nodes, ...
             'FontName', 'Times New Roman', 'FontSize', 24);
    set(gca, 'XTickLabel', repmat({''}, 1, total_nodes));
    set(gca, 'YTickLabel', repmat({''}, 1, total_nodes));
    if any(idx_labels)
        xticks_labels = find(idx_labels);
        yticks_labels = find(idx_labels);
        set(gca, 'XTick', xticks_labels, 'YTick', yticks_labels);
        set(gca, 'XTickLabel', num2str(xticks_labels'));
        set(gca, 'YTickLabel', num2str(yticks_labels'));
    end

end