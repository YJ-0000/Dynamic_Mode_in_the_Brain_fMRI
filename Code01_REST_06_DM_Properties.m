clear; clc; close all;
%%
load DMs/DM_cortical_subcortical_SVD_noROInorm
load results/RSFC_standard_grad

%%
idx_exclude = abs(angle(lambda)) < 1e-10;
lambda(idx_exclude) = [];
Phi_sorted(:,idx_exclude) = [];

Phi_sorted(:,[11,12,17,18]) = Phi_sorted(:,[17,18,11,12]);

N = 360;
label_idx_list = 1:N;

%% FC
N = 360;

for n_dm = 1:6 %1:5
        
    Phi_DM_ang = angle(Phi_sorted(1:N,2*n_dm));
    
    Phi_DM_ang_mat = cos(Phi_DM_ang - Phi_DM_ang');
    
    
    index_tril = tril(ones(size(Phi_DM_ang_mat)),-1);
    Phi_DM_ang_vector = Phi_DM_ang_mat(index_tril ~= 0);
    group_corr_mat_vector = group_corr_mat_original(index_tril ~= 0);
    
    
    figure; scatter(Phi_DM_ang_vector,group_corr_mat_vector);
    
    [r,p] = corrcoef(Phi_DM_ang_vector,group_corr_mat_vector);
    
    disp(r(2));
        
    figure; 
    subplot(1,2,1); heatmap(Phi_DM_ang_mat(1:100,1:100));
    subplot(1,2,2); heatmap(group_corr_mat_original(1:100,1:100));
end

%% lag projection

N = 360;
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

for n_dm = 1:2 %1:5
    Phi_DM_mag = abs(Phi_sorted(1:N,2*n_dm));
        
    Phi_DM_ang = angle(Phi_sorted(1:N,2*n_dm));
    
    Phi_DM_ang_mat = Phi_DM_ang - Phi_DM_ang';
    
    lag_DM = mean(Phi_DM_ang_mat')';
    
    grad1_data = zeros(size(labels.cdata));
    for n_reg = 1:N
        grad1_data(labels.cdata==label_idx_list(n_reg)) = lag_DM(n_reg);
    end
    lag1 = cifti_struct_create_from_template(labels, grad1_data, 'dscalar');
    cifti_write(lag1, ['results/DM',num2str(n_dm),'_lag.dscalar.nii']);
    
    lag_DM_list(:,n_dm) = lag_DM;
end

%% RSFC grad correlation time course
frame_dt = 0.5;
frame_num = 200;
TRtarget = 1.5;
N = 360;

load('results/FPN_vector.mat');
load('results/SMLV_vector.mat');
load('results/DMN_vector.mat');
load('results/Salience_vector.mat');
load('results/DAN_vector.mat');
load('results/VAN_vector.mat');
load('results/LAN_vector.mat');
load('results/Lateralized_index.mat');

for n_dm = 1:6 %1:5
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
    
    Phi_DM_mag = abs(Phi_sorted(1:N,2*n_dm));
    
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    RSFC_grad_temp = RSFC_grad(1:N,:);
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        
        r = corrcoef(eigenstate,RSFC_grad_temp(:,1));
        corr_time_series_PrincipalGradient_DMN(frame+1) = -r(2);
        
        r = corrcoef(eigenstate,RSFC_grad_temp(:,2));
        corr_time_series_SecondGradient(frame+1) = -r(2);
        
        r = corrcoef(eigenstate,DMN_vector);
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
%     legend('Principal Grad','Second Grad','DMN','FPN','salience','SMLV','D attention','V attention','Lateralized index');
    yline(0.7,'k--');


    legend('Principal Grad','Second Grad','DMN','FPN','salience','Language','D attention','V attention','Lateralized index','Location','southeast');
%     legend('Principal Grad','Lateralized index','Location','southeast');
%     legend('Principal Grad','FPN','Salience','threshold');
    xlabel('time (t)');
    ylabel('correlation');
end

%% RSFC grad correlation time course (check)
frame_dt = 0.5;
frame_num = 200;
TRtarget = 1.5;
N = 360;
thres = 0.65;

load('results/FPN_vector.mat');
load('results/SMLV_vector.mat');
load('results/DMN_vector.mat');
load('results/Salience_vector.mat');
load('results/DAN_vector.mat');
load('results/VAN_vector.mat');
load('results/LAN_vector.mat');
load('results/Lateralized_index.mat');

idx_exclude = abs(angle(lambda)) < 1e-10;
lambda(idx_exclude) = [];
Phi_sorted(:,idx_exclude) = [];

for n_dm = 1:(size(Phi_sorted,2)/2) %1:5
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
    
    Phi_DM_mag = abs(Phi_sorted(1:N,2*n_dm));
    
    DM_conjugate1_num = 2*(n_dm-1)+1;
    DM_conjugate2_num = 2*n_dm;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    RSFC_grad_temp = RSFC_grad(1:N,:);
    
    for frame = 0:1:frame_num
        eigenstate = zeros(N,1);
        for roi=1:N 
            eigenstate(roi)= real((lambda_conjugate1^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        
        r = corrcoef(eigenstate,RSFC_grad_temp(:,1));
        corr_time_series_PrincipalGradient_DMN(frame+1) = -r(2);
        
        r = corrcoef(eigenstate,RSFC_grad_temp(:,2));
        corr_time_series_SecondGradient(frame+1) = -r(2);
        
        r = corrcoef(eigenstate,DMN_vector);
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
    end
    
    if max(corr_time_series_PrincipalGradient_DMN)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Grad 1'])
    end
    if max(corr_time_series_SecondGradient)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Grad 2'])
    end
    if max(corr_time_series_DMN)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Default mode (PCC seed)'])
    end
    if max(corr_time_series_D_attention)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Dorsal attention'])
    end
    if max(corr_time_series_V_attention)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Ventral attention'])
    end
    if max(corr_time_series_SMLV)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Sensorimotor seed'])
    end
    if max(corr_time_series_executive)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Executive'])
    end
    if max(corr_time_series_salience)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Salience'])
    end
    if max(corr_time_series_LAN)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Language'])
    end
    if max(corr_time_series_MLI)>thres
        disp(['DM number #',num2str(n_dm,'%03d'), ': Lateralized'])
    end
end
