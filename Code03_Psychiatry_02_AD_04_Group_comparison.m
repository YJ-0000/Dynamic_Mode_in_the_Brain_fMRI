clear; clc;

load DMs/DM_ADNI_cortical_subcortical_noROInorm_indiv_10

%% stack D
abs_all_D = [abs(D_list{1}(2:2:end,:)),...
            abs(D_list{2}(2:2:end,:)),...
            abs(D_list{3}(2:2:end,:)),...
            abs(D_list{4}(2:2:end,:)),...
            abs(D_list{5}(2:2:end,:)),...
            abs(D_list{6}(2:2:end,:))];
angle_all_D = [angle(D_list{1}(2:2:end,:)),...
            angle(D_list{2}(2:2:end,:)),...
            angle(D_list{3}(2:2:end,:)),...
            angle(D_list{4}(2:2:end,:)),...
            angle(D_list{5}(2:2:end,:)),...
            angle(D_list{6}(2:2:end,:))];
        
abs_all_D = abs_all_D';
angle_all_D = angle_all_D';
%% make regressors
num_all_sub = size(D_list{1},2)+size(D_list{2},2) + size(D_list{3},2)+size(D_list{4},2) + size(D_list{5},2)+size(D_list{6},2);
is_AD = [ones(size(D_list{1},2)+size(D_list{2},2),1);
        zeros(size(D_list{3},2)+size(D_list{4},2),1);
        zeros(size(D_list{5},2)+size(D_list{6},2),1)];
is_MCI = [zeros(size(D_list{1},2)+size(D_list{2},2),1);
        zeros(size(D_list{3},2)+size(D_list{4},2),1);
        ones(size(D_list{5},2)+size(D_list{6},2),1)];
is_CN = [zeros(size(D_list{1},2)+size(D_list{2},2),1);
        ones(size(D_list{3},2)+size(D_list{4},2),1);
        zeros(size(D_list{5},2)+size(D_list{6},2),1)];

%%
all_func_seq = [AD_sub_info{1}.func_seq;
                AD_sub_info{2}.func_seq;
                AD_sub_info{3}.func_seq;
                AD_sub_info{4}.func_seq;
                AD_sub_info{5}.func_seq;
                AD_sub_info{6}.func_seq];

all_func_manu = [AD_sub_info{1}.func_manu;
                AD_sub_info{2}.func_manu;
                AD_sub_info{3}.func_manu;
                AD_sub_info{4}.func_manu;
                AD_sub_info{5}.func_manu;
                AD_sub_info{6}.func_manu];
            
func_seq_types = unique(all_func_seq);
func_manu_types = unique(all_func_manu);

%%
func_seq_regressor = zeros(num_all_sub,length(func_seq_types));
func_manu_regressor = zeros(num_all_sub,length(func_manu_types));
n_accum_sub = 0;
for ndir = 1:length(sub_list)
    temp_sub_list = sub_list{ndir};
    sub_info = AD_sub_info{ndir};
    for nsub = 1:length(temp_sub_list)
        sub_name = temp_sub_list{nsub};
%         fprintf(['sub >> ',sub_name,'...']);
        
        n_accum_sub = n_accum_sub  + 1;
        sub_idx = find(strcmp(sub_info.subject_id',sub_name));
        
        seq_name = sub_info(sub_idx,'func_seq').Variables;
        seq_idx = find(strcmp(func_seq_types,seq_name));
        func_seq_regressor(n_accum_sub,seq_idx) = 1;
        
        manu_name = sub_info(sub_idx,'func_manu').Variables;
        manu_idx = find(strcmp(func_manu_types,manu_name));
        func_manu_regressor(n_accum_sub,manu_idx) = 1;
        
        
%         fprintf('Done.\n');
    end
end
%% T-test
disp('########### T-test ###########')
% Define explanatory variables
X = [is_AD,is_MCI,is_CN];
% Define predictor variables (D)
y_abs = abs_all_D;
y_angle = angle_all_D;

% Define confounding regressors (func_seq_regressor and func_manu_regressor)
confound_regressors = [func_seq_regressor, func_manu_regressor];

% Define contrast vectors
c_AD_minus_CN = zeros(size(X,2)+size(confound_regressors,2),1);
c_AD_minus_CN(1) = 1; c_AD_minus_CN(3) = -1;
c_MCI_minus_CN = zeros(size(X,2)+size(confound_regressors,2),1);
c_MCI_minus_CN(2) = 1; c_MCI_minus_CN(3) = -1;
c_AD_minus_MCI = zeros(size(X,2)+size(confound_regressors,2),1);
c_AD_minus_MCI(1) = 1; c_AD_minus_MCI(2) = -1;

%%% test magnitude
for n_dm = 1:5
    % Calculate T-value 
    [T_value_AD_minus_CN,p_AD_minus_CN] = t_test_GLM([X,confound_regressors],y_abs(:,n_dm),c_AD_minus_CN);
    [T_value_MCI_minus_CN,p_MCI_minus_CN] = t_test_GLM([X,confound_regressors],y_abs(:,n_dm),c_MCI_minus_CN);
    [T_value_AD_minus_MCI,p_AD_minus_MCI] = t_test_GLM([X,confound_regressors],y_abs(:,n_dm),c_AD_minus_MCI);

    % Display T-values
    fprintf('DM #%d (magnitude) -- T-value for AD-CN: T=%.4f, p=%.4f\n', n_dm, T_value_AD_minus_CN,p_AD_minus_CN);
    fprintf('DM #%d (magnitude) -- T-value for MCI-CN: T=%.4f, p=%.4f\n', n_dm, T_value_MCI_minus_CN,p_MCI_minus_CN);
    fprintf('DM #%d (magnitude) -- T-value for AD-MCI: T=%.4f, p=%.4f\n', n_dm, T_value_AD_minus_MCI,p_AD_minus_MCI);
    
    % Calculate T-value 
    [T_value_AD_minus_CN,p_AD_minus_CN] = t_test_GLM([X,confound_regressors],y_angle(:,n_dm),c_AD_minus_CN);
    [T_value_MCI_minus_CN,p_MCI_minus_CN] = t_test_GLM([X,confound_regressors],y_angle(:,n_dm),c_MCI_minus_CN);
    [T_value_AD_minus_MCI,p_AD_minus_MCI] = t_test_GLM([X,confound_regressors],y_angle(:,n_dm),c_AD_minus_MCI);

    % Display T-values
    fprintf('DM #%d (phase) -- T-value for AD-CN: T=%.4f, p=%.4f\n', n_dm, T_value_AD_minus_CN,p_AD_minus_CN);
    fprintf('DM #%d (phase) -- T-value for MCI-CN: T=%.4f, p=%.4f\n', n_dm, T_value_MCI_minus_CN,p_MCI_minus_CN);
    fprintf('DM #%d (phase) -- T-value for AD-MCI: T=%.4f, p=%.4f\n', n_dm, T_value_AD_minus_MCI,p_AD_minus_MCI);
end
%% F-test
disp('########### F-test ###########')
%%% test magnitude
for n_dm = 1:5
    [F,p_val] = F_test_GLM([X,confound_regressors],y_abs(:,n_dm),[c_AD_minus_MCI,c_MCI_minus_CN]);
    fprintf('DM #%d (magnitude) -- F-value for AD vs MCI vs CN: F=%.4f, p=%.4f\n', n_dm, F,p_val);
    [F,p_val] = F_test_GLM([X,confound_regressors],y_angle(:,n_dm),[c_AD_minus_MCI,c_MCI_minus_CN]);
    fprintf('DM #%d (phase) -- F-value for AD vs MCI vs CN: F=%.4f, p=%.4f\n', n_dm, F,p_val);
end

%% statistical function
function [t,p_val] = t_test_GLM(X,y,c,d)
    if nargin < 4
        d = 0;
    end
    
    J = size(X,1);
    p = rank(X);
    
    beta = pinv(X) * y;
    
    resid_e = y - X * beta;
    
    sigma_squared = (resid_e' * resid_e)/(J-p);
    
    t = (c'*beta - d)/sqrt(sigma_squared * c' * pinv(X'*X) * c); %#ok<MINV>
    p_val = 2 * (1 - tcdf(abs(t), J-p));
end

function [F,p_val] = F_test_GLM(X,y,C)
    J = size(X,1);
    p = rank(X);
    
    C0 = eye(size(X,2)) - C*pinv(C);
    
    X0 = X * C0;
    
    p2 = rank(X0);
    
    p1 = p - p2;
    
    
    R0 = eye(J) - X0 * pinv(X0);
    R = eye(J) - X * pinv(X);
    
    M = R0 - R;
    
    F = ((J-p)/p1) * (y' * M * y) / (y' * R * y);
    p_val = 1 - fcdf(F, p1, J-p);
end