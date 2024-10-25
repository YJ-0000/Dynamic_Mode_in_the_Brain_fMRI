clear; clc;

load DMs/DM_MDD_cortical_subcortical_noROInorm_indiv_10

%% stack D
abs_all_D = [abs(D_list{1}(2:2:end,:)),...
            abs(D_list{2}(2:2:end,:))];
angle_all_D = [angle(D_list{1}(2:2:end,:)),...
            angle(D_list{2}(2:2:end,:))];
        
abs_all_D = abs_all_D';
angle_all_D = angle_all_D';
%% make regressors
num_all_sub = size(D_list{1},2)+size(D_list{2},2);
is_MDD = [zeros(size(D_list{1},2),1);
        ones(size(D_list{2},2),1)];
is_HC = [ones(size(D_list{1},2),1);
        zeros(size(D_list{2},2),1)];

%%

site_list = {'SWA','HUH','HRC','HKH','COI','KUT','KTT','UTO','ATT','ATV','CIN','OSU','NKN'};

%%
site_regressor = zeros(num_all_sub,length(site_list));
age_regressor = zeros(num_all_sub,1);
sex_regressor = zeros(num_all_sub,1);
n_accum_sub = 0;
for ndir = 1:length(sub_list)
    temp_sub_list = sub_list{ndir};
    for nsub = 1:length(temp_sub_list)
        sub_name = temp_sub_list{nsub};
%         fprintf(['sub >> ',sub_name,'...']);
        
        n_accum_sub = n_accum_sub  + 1;
        sub_idx = find(strcmp(sub_info.participant_id,sub_name));
        assert(~isempty(sub_idx));
        
        site_name = sub_info(sub_idx,'site').Variables;
        site_idx = find(strcmp(site_list,site_name));
        
        site_regressor(n_accum_sub,site_idx) = 1;
        
        age_regressor(n_accum_sub) = sub_info(sub_idx,'age').Variables;
        sex_regressor(n_accum_sub) = sub_info(sub_idx,'sex').Variables;
        
        
%         fprintf('Done.\n');
    end
end
site_regressor(:,sum(site_regressor)==0) = [];
%% T-test
disp('########### T-test ###########')
% Define explanatory variables
X = [is_MDD,is_HC];
% Define predictor variables (D)
y_abs = abs_all_D;
y_angle = angle_all_D;

% Define confounding regressors (func_seq_regressor and func_manu_regressor)
confound_regressors = [site_regressor];

% Define contrast vectors
c_MDD_minus_HC = zeros(size(X,2)+size(confound_regressors,2),1);
c_MDD_minus_HC(1) = 1; c_MDD_minus_HC(2) = -1;

%%% test magnitude
for n_dm = 1:5
    % Calculate T-value 
    [T_value_MDD_minus_HC,p_MDD_minus_HC] = t_test_GLM([X,confound_regressors],y_abs(:,n_dm),c_MDD_minus_HC);

    % Display T-values
    fprintf('DM #%d (magnitude) -- T-value for MDD-HC: T=%.4f, p=%.4f\n', n_dm, T_value_MDD_minus_HC,p_MDD_minus_HC);
    
    % Calculate T-value 
    [T_value_MDD_minus_HC,p_MDD_minus_HC] = t_test_GLM([X,confound_regressors],y_angle(:,n_dm),c_MDD_minus_HC);

    % Display T-values
    fprintf('DM #%d (phase) -- T-value for MDD-HC: T=%.4f, p=%.4f\n', n_dm, T_value_MDD_minus_HC,p_MDD_minus_HC);
end
%% F-test
disp('########### F-test ###########')
%%% test magnitude
% for n_dm = 1:5
%     [F,p_val] = F_test_GLM([X,confound_regressors],y_abs(:,n_dm),[c_AD_minus_MCI,c_MCI_minus_CN]);
%     fprintf('DM #%d (magnitude) -- F-value for AD vs MCI vs CN: T=%.4f, p=%.4f\n', n_dm, F,p_val);
%     [F,p_val] = F_test_GLM([X,confound_regressors],y_angle(:,n_dm),[c_AD_minus_MCI,c_MCI_minus_CN]);
%     fprintf('DM #%d (phase) -- F-value for AD vs MCI vs CN: T=%.4f, p=%.4f\n', n_dm, F,p_val);
% end

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