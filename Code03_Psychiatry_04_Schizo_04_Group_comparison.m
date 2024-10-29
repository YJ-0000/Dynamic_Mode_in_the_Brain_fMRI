clear; clc;
current_path = pwd;
load DMs/DM_Schizo_cortical_subcortical_noROInorm_indiv_10

%% stack D
abs_all_D = abs(D(2:2:end,:));
angle_all_D = angle(D(2:2:end,:));
abs_all_B = (B(1:2:end,:));
        
abs_all_D = abs_all_D';
angle_all_D = angle_all_D';
abs_all_B = abs_all_B';
%% make regressors

data_path = 'G:\COBRE_preprocessed\4197885';
cd(data_path);
sub_info_table = readtable('phenotypic_data_sorted.txt','FileType', 'text', 'Delimiter', '\t');

age_regressor = sub_info_table.CurrentAge;
sex_regressor = double(strcmp(sub_info_table.Gender,'Female'));
handedness_regressor = double(strcmp(sub_info_table.Handedness,'Right'));
FD_regressor = sub_info_table.FD;
FramesOK_regressor = sub_info_table.FramesOK;
FDScrubbed_regressor = sub_info_table.FDScrubbed;

is_patient = double(strcmp(sub_info_table.SubjectType,'Patient'));
is_control = double(strcmp(sub_info_table.SubjectType,'Control'));

cd(current_path);
%% T-test
disp('########### T-test ###########')
% Define explanatory variables
X = [is_patient,is_control];
% Define predictor variables (D)
y_abs = abs_all_D;
y_angle = angle_all_D;
y_B = abs_all_B;

% Define confounding regressors (func_seq_regressor and func_manu_regressor)
confound_regressors = [age_regressor,sex_regressor,handedness_regressor,FD_regressor,FramesOK_regressor,FDScrubbed_regressor];

% Define contrast vectors
c_Schizo_minus_HC = zeros(size(X,2)+size(confound_regressors,2),1);
c_Schizo_minus_HC(1) = 1; c_Schizo_minus_HC(2) = -1;

%%% test magnitude
for n_dm = 1:5
    % Calculate T-value 
    [T_value_Schizo_minus_HC,p_Schizo_minus_HC] = t_test_GLM([X,confound_regressors],y_abs(:,n_dm),c_Schizo_minus_HC);

    % Display T-values
    fprintf('DM #%d (magnitude) -- T-value for Schizo-HC: T=%.4f, p=%.4f\n', n_dm, T_value_Schizo_minus_HC,p_Schizo_minus_HC);
    
    % Calculate T-value 
    [T_value_Schizo_minus_HC,p_Schizo_minus_HC] = t_test_GLM([X,confound_regressors],y_angle(:,n_dm),c_Schizo_minus_HC);

    % Display T-values
    fprintf('DM #%d (phase) -- T-value for Schizo-HC: T=%.4f, p=%.4f\n', n_dm, T_value_Schizo_minus_HC,p_Schizo_minus_HC);
    
    % Calculate T-value 
    [T_value_Schizo_minus_HC,p_Schizo_minus_HC] = t_test_GLM([X,confound_regressors],y_B(:,n_dm),c_Schizo_minus_HC);

    % Display T-values
    fprintf('DM #%d (Intensity) -- T-value for Schizo-HC: T=%.4f, p=%.4f\n', n_dm, T_value_Schizo_minus_HC,p_Schizo_minus_HC);
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