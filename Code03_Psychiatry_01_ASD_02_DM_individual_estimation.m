% clear;clc;

%%
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm
load results/ABIDE_timeseires_cortical_subcortical_CPAC_denoised_filtered

siteList={'CALTECH','CMU','KKI','MAXMUN','NYU','OLIN','OHSU','SDSU','SBL','STANFORD','TRINITY','UCLA','LEUVEN','UM','PITT','USM','YALE'};
TRlist=[2000 2000 2500 3000 2000 1500 2500 2000 2200 2000 2000 3000 1656 2000 1500 2000 2000];
%% individual fitting
TRtarget = 1.5;
thres = 1;

num_DMs = 10;
U=Phi_sorted(:,1:num_DMs);
count_roi = zeros(size(U,1),1);
for n=1:length(time_series_denoised_filtered)
    if isempty(time_series_denoised_filtered{n}) 
        continue;
    end
    y = time_series_denoised_filtered{n};
    y(roi_exclude,:) = [];
    var_y = var(y,0,2);
    count_roi = count_roi + (var_y<0.0001);
end
U(count_roi>thres,:)=[];
V=pinv(U);
residual_matrix = eye(size(U,1)) - U(:,1:num_DMs)*V(1:num_DMs,:);
% V(:,count_roi>thres) = [];
D=zeros(num_DMs+1,length(time_series_denoised_filtered));
B=zeros(num_DMs,length(time_series_denoised_filtered));
idx_exclude = false(1,length(time_series_denoised_filtered));
for n=1:length(time_series_denoised_filtered)
    if isempty(time_series_denoised_filtered{n}) 
        idx_exclude(n) = true;
        continue;
    end
    
    disp(['start: sub#' num2str(n)]);
    Z_temp = zeros(num_DMs+1,num_DMs+1);
    W_temp = zeros(num_DMs+1,1);
    
    tic
    
    index = 1;
    while ~contains(upper(image_file_list(n).name),siteList{index})
        index = index + 1;
    end
    
    y = time_series_denoised_filtered{n};
%     y = y(1:360,:);
    y(roi_exclude,:) = [];
    y(count_roi>thres,:)=[];
    t = (1:size(y,2)) * (TRlist(index)/1000);
    t_fine = TRtarget:TRtarget:t(end);
    if TRlist(index)/1000 ~= TRtarget
        pp = spline(t, y);          % Compute piecewise polynomial (B-spline) representation
        y_fine = ppval(pp, t_fine); % Evaluate the piecewise polynomial at the finer time points
    else
        y_fine = y;
    end
    
    X_temp = y_fine(:,2:end);
    Y_temp = y_fine(:,1:end-1);
    resid_Y_temp = residual_matrix * Y_temp;
    VY = V(1:num_DMs,:) * Y_temp;
    Z_temp(1,1) = sum(dot(resid_Y_temp, resid_Y_temp, 1));
    for j=1:num_DMs
        Z_temp(1,j+1) = sum(dot(U(:,j)'*resid_Y_temp, VY(j,:), 1));
    end
    Z_temp(2:end,1) = Z_temp(1,2:end)';
    for i=1:num_DMs
        for j=1:num_DMs  
            Z_temp(i+1,j+1) = (U(:,i)'*U(:,j))*sum(dot(VY(i,:), VY(j,:), 1));
        end
    end
    W_temp(1) = sum(dot(resid_Y_temp,X_temp,1));
    for k=1:num_DMs
        W_temp(k+1) = sum(dot(U(:,k)*VY(k,:),X_temp,1));
    end
    toc

    D(:,n) = Z_temp\W_temp;
    B(:,n) = mean(abs(VY),2);
%     D(:,n) = C_temp(2:end,2:end)\B_temp(2:end);
    disp(['end: sub#' num2str(n)]);
end
%%
D(:,idx_exclude) = [];
B(:,idx_exclude) = [];
image_file_list(idx_exclude) = [];

save DMs/DM_ABIDE_cortical_subcortical_noROInorm_indiv_10 D B image_file_list Phi_sorted lambda idx_exclude