clear; clc;
current_path = pwd;
close all;
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered_meta.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered.mat');


%%

load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm

n_time = 1200;

t_sample = 0.72;
% TRtarget = 0.72;
TRtarget = 1.5;

t = (1:n_time) * (t_sample);
t_fine = TRtarget:TRtarget:t(end);

num_DMs = 18;
U=Phi_sorted(:,1:num_DMs);
V=pinv(U);
D1=zeros(num_DMs+1,size(time_series_denoised_filtered,1));
D2=zeros(num_DMs+1,size(time_series_denoised_filtered,1));
idx_exclude = false(1,size(time_series_denoised_filtered,1));
for nsub = 1:size(time_series_denoised_filtered,1)
    tau_temp = 0;
    disp(nsub);
    X1 = zeros(N-sum(roi_exclude),2*(length(t_fine)-1));
    Y1 = X1;
    X2 = X1;
    Y2 = X1;
    n_accum = 0;
    time_series_denoised_filtered_temp = time_series_denoised_filtered(nsub,[1,2,3,4]);
    for nses = 1:4
        if ~isempty(time_series_denoised_filtered_temp{1,nses})
            y = time_series_denoised_filtered_temp{1,nses};
            y(roi_exclude,:) = [];
            if t_sample ~= TRtarget
                pp = spline(t, y);          % Compute piecewise polynomial (B-spline) representation
                y_fine = ppval(pp, t_fine); % Evaluate the piecewise polynomial at the finer time points
            else
                y_fine = y;
            end
            if nses <= 2
                X1(:,(nses-1)*(length(t_fine)-1)+1:nses*(length(t_fine)-1)) = y_fine(:,2:end);
                Y1(:,(nses-1)*(length(t_fine)-1)+1:nses*(length(t_fine)-1)) = y_fine(:,1:end-1);
            else
                X2(:,(nses-3)*(length(t_fine)-1)+1:(nses-2)*(length(t_fine)-1)) = y_fine(:,2:end);
                Y2(:,(nses-3)*(length(t_fine)-1)+1:(nses-2)*(length(t_fine)-1)) = y_fine(:,1:end-1);
            end
            n_accum = n_accum + 1;
        else
            idx_exclude(nsub) = true;
            break
        end
    end
    
    tic
    if n_accum == 4
        for n_part = 1:2
            C_temp = zeros(num_DMs+1,num_DMs+1);
            B_temp = zeros(num_DMs,1);

            if n_part == 1
                X_temp = X1;
                Y_temp = Y1;
            else
                X_temp = X2;
                Y_temp = Y2;
            end
            VY = V(1:num_DMs,:) * Y_temp;
            C_temp(1,1) = sum(dot(Y_temp, Y_temp, 1));
            for j=1:num_DMs
                C_temp(1,j+1) = sum(dot(U(:,j)'*Y_temp, VY(j,:), 1));
            end
            C_temp(2:end,1) = C_temp(1,2:end)';
            for i=1:num_DMs
                for j=1:num_DMs  
                    C_temp(i+1,j+1) = (U(:,i)'*U(:,j))*sum(dot(VY(i,:), VY(j,:), 1));
                end
            end
            B_temp(1) = sum(dot(Y_temp,X_temp,1));
            for k=1:num_DMs
                B_temp(k+1) = sum(dot(U(:,k)*VY(k,:),X_temp,1));
            end
            if n_part == 1
                D1(:,nsub) = C_temp\B_temp;
            else
                D2(:,nsub) = C_temp\B_temp;
            end
        end
    end
    toc
    
end

D1(:,idx_exclude) = [];
D2(:,idx_exclude) = [];
sub_ids(idx_exclude) = [];

save DMs/DMs_REST_test_retest_indiv_18 D1 D2 sub_ids

%%
rest_self_corr = zeros(1,num_DMs);
rest_self_sig = zeros(1,num_DMs);
% rest_self_ICC = zeros(1,10);
for n_dm = 1:(num_DMs/2)+1
    temp_D1 = D1(2*(n_dm-1)+1,:);
    temp_D2 = D2(2*(n_dm-1)+1,:);
    
    [r,p] = corrcoef(abs(temp_D1),abs(temp_D2));
    rest_self_corr(2*(n_dm-1)+1) = r(2);
    rest_self_sig(2*(n_dm-1)+1) = p(2);
    [r,p] = corrcoef(angle(temp_D1),angle(temp_D2));
    rest_self_corr(2*(n_dm)) = r(2);
    rest_self_sig(2*(n_dm)) = p(2);
    
%     rest_self_ICC(2*(n_dm-1)+1) = computeICC(abs(temp_D1),abs(temp_D2));
%     rest_self_ICC(2*(n_dm)) = computeICC(angle(temp_D1),angle(temp_D2));
end
figure; 
heatmap(rest_self_corr);
