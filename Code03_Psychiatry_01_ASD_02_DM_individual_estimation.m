clear;clc;

%%
load DMs/DM_cortical_subcortical_SVD_noROInorm
load results/ABIDE_timeseires_cortical_subcortical_CPAC_denoised_filtered

siteList={'CALTECH','CMU','KKI','MAXMUN','NYU','OLIN','OHSU','SDSU','SBL','STANFORD','TRINITY','UCLA','LEUVEN','UM','PITT','USM','YALE'};
TRlist=[2000 2000 2500 3000 2000 1500 2500 2000 2200 2000 2000 3000 1656 2000 1500 2000 2000];
%% individual fitting
TRtarget = 1.5;

idx_exclude = abs(angle(lambda)) < 1e-10;
lambda(idx_exclude) = [];
Phi_sorted(:,idx_exclude) = [];

max_number_eigenstates = 12;
U=Phi_sorted;
U=U(1:360,:);
% U(:,[9,10,11,12,17,18]) = U(:,[11,12,17,18,9,10]);
U(:,[11,12,17,18]) = U(:,[17,18,11,12]);
V=pinv(U);
D=zeros(max_number_eigenstates,length(time_series_denoised_filtered));
idx_exclude = false(1,length(time_series_denoised_filtered));
for n=1:length(time_series_denoised_filtered)
    if isempty(time_series_denoised_filtered{n}) 
        idx_exclude(n) = true;
        continue;
    end
    
    disp(['start: sub#' num2str(n)]);
    C_temp = zeros(max_number_eigenstates,max_number_eigenstates);
    B_temp = zeros(max_number_eigenstates,1);
    
    tic
    
    index = 1;
    while ~contains(upper(image_file_list(n).name),siteList{index})
        index = index + 1;
    end
    
    y = time_series_denoised_filtered{n};
    y = y(1:360,:);
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
    VY = V(1:max_number_eigenstates,:) * Y_temp;
    for i=1:max_number_eigenstates
        for j=1:max_number_eigenstates  
            C_temp(i,j) = (U(:,i)'*U(:,j))*sum(dot(VY(i,:), VY(j,:), 1));
        end
    end
    for k=1:max_number_eigenstates
        B_temp(k) = sum(dot(U(:,k)*VY(k,:),X_temp,1));
    end
    toc

    D(:,n) = C_temp\B_temp;
    disp(['end: sub#' num2str(n)]);
end
%%
D(:,idx_exclude) = [];
image_file_list(idx_exclude) = [];

save DMs/DM_ABIDE_cortical_subcortical_SVD_noROInorm_indiv_12 D image_file_list Phi_sorted lambda idx_exclude