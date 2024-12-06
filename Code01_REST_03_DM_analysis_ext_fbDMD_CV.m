clear; clc;
current_path = pwd;
% conn;
close all;
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered_meta.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered.mat');
%%

n_time = 1200;
i_num = 0;
for ii = 1:(4*size(time_series_denoised_filtered,1))
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if isempty(time_series_denoised_filtered{nsub,nses})
        continue
    end
    i_num = i_num + 1;
end
disp(i_num)

t_sample = 0.72;
% TRtarget = 0.72;
TRtarget = 1.5;

t = (1:n_time) * (t_sample);
t_fine = TRtarget:TRtarget:t(end);

X = zeros(N,i_num * (length(t_fine)-1));
Y = X;
data_index_sub = zeros(1,i_num * (length(t_fine)-1),'uint16');

tau_index_sub = zeros(i_num,1);

num_sub = size(time_series_denoised_filtered,1);
tau=zeros(num_sub,1);

i_num = 0;
for nsub = 1:size(time_series_denoised_filtered,1)
    tau_temp = 0;
    disp(nsub);
    for nses = 1:4
        if ~isempty(time_series_denoised_filtered{nsub,nses})
            y = time_series_denoised_filtered{nsub,nses};
            if t_sample ~= TRtarget
                pp = spline(t, y);          % Compute piecewise polynomial (B-spline) representation
                y_fine = ppval(pp, t_fine); % Evaluate the piecewise polynomial at the finer time points
            else
                y_fine = y;
            end
            
            if isnan(sum(y_fine,'all'))
                warning('There is NAN!!')
            end
            
            i_num = i_num + 1;
            tau_index_sub(i_num) = nsub;

            tau_temp = tau_temp + (size(y_fine,2)-1);
            X(:,(i_num-1)*(length(t_fine)-1)+1:i_num*(length(t_fine)-1)) = y_fine(:,2:end);
            Y(:,(i_num-1)*(length(t_fine)-1)+1:i_num*(length(t_fine)-1)) = y_fine(:,1:end-1);
            data_index_sub((i_num-1)*(length(t_fine)-1)+1:i_num*(length(t_fine)-1)) = nsub;
        end
    end
    tau(nsub) = tau_temp;
end
X(:,i_num*(length(t_fine)-1)+1:end) = [];
Y(:,i_num*(length(t_fine)-1)+1:end) = [];
tau_index_sub(i_num+1:end) = [];

disp('*** total number of time points ***');
disp([i_num*(length(t_fine)-1), size(X,2)]);

var_X = var(X,0,2);
roi_exclude = var_X < 0.001;
X(roi_exclude,:) = [];
Y(roi_exclude,:) = [];

clear time_series_denoised_filtered

%% Cross-validation parametters
disp('*** CV initialize ***');
cv_num = 5;
rng(123);
sub_id_perm = randperm(num_sub);
train_id_list = cell(cv_num,1);
test_id_list = cell(cv_num,1);
% Calculate the number of subjects per fold
fold_size = floor(num_sub / cv_num);

% Loop through each fold to assign training and testing IDs
for n_cv = 1:cv_num
    % Define the start and end indices for the current test fold
    start_idx = (n_cv - 1) * fold_size + 1;
    
    if n_cv < cv_num
        end_idx = n_cv * fold_size;
    else
        % Ensure the last fold includes any remaining subjects
        end_idx = num_sub;
    end
    
    % Extract the test IDs for the current fold
    test_id = sub_id_perm(start_idx:end_idx);
    
    % The training IDs are all other IDs not in the current test fold
    train_id = setdiff(sub_id_perm, test_id);
    
    % Store the IDs in the respective cell arrays
    test_id_list{n_cv} = sort(test_id);
    train_id_list{n_cv} = train_id;
end

%% Cross-validation 
disp('*** CV DM estimation started! ***');
lambda_fold = cell(cv_num,1);
Phi_fold = cell(cv_num,1);
Phi_all_fold = cell(cv_num,1);
lambda_all_fold = cell(cv_num,1);
for n_cv = 1:cv_num
    disp(['=== CV DM estimation fold #',num2str(n_cv),'===']);
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    is_train = ismember(data_index_sub, current_train_subjects);
    is_test = ismember(data_index_sub, current_test_subjects);
    disp(['Data points in X_train: ', num2str(sum(is_train))]);
    
    clear X_train Y_train
    Y_train = Y(:,is_train);
    X_train = X(:,is_train);
    
    %%% fbDMD
    disp(['> CV fold #',num2str(n_cv),': Extended fbDMD']);
    tic
    A1 = X_train*Y_train'; A2 = Y_train*Y_train';
    A_f = A1 * pinv(A2);
    B1 = Y_train*X_train'; B2 = X_train*X_train';
    clear X_train Y_train
    A_b = B1 * pinv(B2);
    A = (A_f/A_b)^0.5;
    A = real(A);
    
    [Phi_sorted,D] = eig(A);
    lambda = diag(D);
    idx_exclude = (abs(angle(lambda)) < 2*pi*1.5*0.01) | (abs(angle(lambda)) > 2*pi*1.5*0.1);
    lambda_rest  = lambda(idx_exclude);
    lambda(idx_exclude) = [];
    Phi_rest = Phi_sorted(:,idx_exclude);
    Phi_sorted(:,idx_exclude) = [];
    [lambda,idx_sort] = sort(lambda,'descend');
    Phi_sorted = Phi_sorted(:,idx_sort);
    [lambda_rest,idx_sort] = sort(lambda_rest,'descend');
    Phi_rest_sorted = Phi_rest(:,idx_sort);
    
    toc
    disp(['=== DONE ===']);
    
    lambda_fold{n_cv} = lambda;
    Phi_fold{n_cv} = Phi_sorted;
    Phi_all_fold{n_cv} = [Phi_sorted,Phi_rest_sorted];
    lambda_all_fold{n_cv} = [lambda;lambda_rest];
end

%% save CV results
save DMs/DM_cortical_subcortical_noROInorm_CV Phi_fold lambda_fold roi_exclude

%% Testing DM consistency
max_DMs = 24;


rng(321); % for reproducibility

Phi_raw_all = [];
Phi_abs_all = [];
Phi_angle_all = [];
for n_cv = 1:cv_num
    Phi = Phi_fold{n_cv};
    Phi_raw_all = [Phi_raw_all,(Phi(:,1:2:max_DMs))]; %#ok<AGROW>
    Phi_abs_all = [Phi_abs_all,abs(Phi(:,1:2:max_DMs))]; %#ok<AGROW>
    Phi_angle_all = [Phi_angle_all,angle(Phi(:,1:2:max_DMs))]; %#ok<AGROW>
end
[cluster_idx, C] = kmeans(Phi_abs_all', max_DMs/2, 'Replicates', 10);
for n_cv = 1:cv_num
    assert(length(unique(cluster_idx((n_cv-1)*max_DMs/2+1:n_cv*max_DMs/2))) == max_DMs/2)
end

consistency_corr_list = cell(max_DMs/2,1);
consistency_CS_list = cell(max_DMs/2,1);
for n_dm = 1:max_DMs/2
    figure;
    iter_n = 0;
    Phi_select = Phi_abs_all(:,cluster_idx==n_dm);
    Phi_raw_select = Phi_raw_all(:,cluster_idx==n_dm);
    corr_mat = zeros(cv_num);
    CS_mat = zeros(cv_num);
    for n_cv1 = 1:cv_num
        for n_cv2 = 1:cv_num
            iter_n = iter_n + 1;
            subplot(cv_num,cv_num,iter_n);
            scatter(Phi_select(:,n_cv1),Phi_select(:,n_cv2))
            
            r = corrcoef(Phi_select(:,n_cv1),Phi_select(:,n_cv2));
            corr_mat(n_cv1,n_cv2) = r(2);
            
            cs = abs(dot(Phi_raw_select(:,n_cv1),Phi_raw_select(:,n_cv2)))/...
                (sqrt(dot(Phi_raw_select(:,n_cv1),Phi_raw_select(:,n_cv1))) * sqrt(dot(Phi_raw_select(:,n_cv2),Phi_raw_select(:,n_cv2))));

            CS_mat(n_cv1,n_cv2) = abs(cs);
        end
    end
    consistency_corr_list{n_dm} = corr_mat;
    consistency_CS_list{n_dm} = CS_mat;
end

DM_labels = {'Principal DM (DM 1)','SN-to-DMN (DM 2)','SN-to-CEN (DM 3)','FV-to-SM (DM 4)','Bi-asym (DM 5)', ...
             'DM 6','DM 7', 'DM 8', 'DM 9', ...
             'DM 10', 'DM 11', 'DM 12'};
            
figure;
count = 0;
for n_dm = [3,10,2,5,7,12,8,9,1,11,6,4]
    count = count + 1;
    subplot(3,4,count); 
    heatmap(consistency_CS_list{n_dm}, 'ColorLimits', [0,1]);
    title(DM_labels{count});
    if rem(count,4) == 1; ylabel('CV fold'); end
    if count > 8; xlabel('CV fold'); end
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
end

% figure;
% count = 0;
% for n_dm = [3,10,2,5,7,12,8,9,1,11,6,4]
%     count = count + 1;
%     subplot(3,4,count); 
%     heatmap(consistency_corr_list{n_dm}, 'ColorLimits', [0,1]);
%     title(DM_labels{count});
%     if rem(count,4) == 1; ylabel('CV fold'); end
%     if count > 8; xlabel('CV fold'); end
%     set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
% end
%% CV prediction
disp('*** CV prediction started! ***');
max_DMs = 24;
loss_fold = zeros(cv_num,max_DMs/2);
R_squared_zero_fold = zeros(cv_num,max_DMs/2);
R_squared_null_fold = zeros(cv_num,max_DMs/2);
R_squared_fold = zeros(cv_num,max_DMs/2);
for n_cv = 1:cv_num
    disp(['=== CV prediction fold #',num2str(n_cv),'===']);
    lambda = lambda_fold{n_cv};
    Phi_sorted = Phi_fold{n_cv};
    Phi_all = Phi_all_fold{n_cv};
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    for num_DMs = 2:2:max_DMs
        tic
        
        %%% Projecting BOLD acitivity to the target sub-space not considering other possible modes
        U=Phi_sorted(:,1:num_DMs);
        V=pinv(U);
        residual_matrix = eye(size(U,1)) - U*V;

        %%% Projecting BOLD acitivity to the target sub-space only considering physically plausible modes
%         U=Phi_sorted;
%         V=pinv(U);
%         residual_matrix = eye(size(U,1)) - U(:,1:num_DMs)*V(1:num_DMs,:);
        
        %%% Projecting BOLD acitivity to the target sub-space considering all group-level modes
%         U=Phi_all;
%         V=pinv(U);
%         residual_matrix = eye(size(U,1)) - U(:,1:num_DMs)*V(1:num_DMs,:);
        
        loss_fold_current = zeros(1,length(current_test_subjects));
        R_squared_fold_current = loss_fold_current;
        R_squared_fold_null_current = loss_fold_current;
        R_squared_fold_zero_current = loss_fold_current;
        session_count = zeros(1,length(current_test_subjects));
        for n_test = 1:length(current_test_subjects)
            current_subject = current_test_subjects(n_test);
            is_test = data_index_sub == current_subject;
%             disp(['Data points in X_train: ', num2str(sum(is_test))]);

            num_session = round(sum(is_test)/(length(t_fine)-1));
            X_test = X(:,is_test);
            Y_test = Y(:,is_test);
            loss_session = zeros(1,num_session);
            R_squared_zero_session = loss_session;
            R_squared_null_session = loss_session;
            R_squared_session = loss_session;
            for n_ses = 1:num_session
                X_test_sess = X_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
                Y_test_sess = Y_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
                C_temp = zeros(num_DMs+1,num_DMs+1);
                B_temp = zeros(num_DMs+1,1);

                X_temp = X_test_sess(:,1:round(length(t_fine)*4/5));
                Y_temp = Y_test_sess(:,1:round(length(t_fine)*4/5));
                resid_Y_temp = residual_matrix * Y_temp;
                X_temp_pred = X_test_sess(:,round(length(t_fine)*4/5)+1:end);
                Y_temp_pred = Y_test_sess(:,round(length(t_fine)*4/5)+1:end);
                VY = V(1:num_DMs,:) * Y_temp;
                C_temp(1,1) = sum(dot(resid_Y_temp, resid_Y_temp, 1));
                for j=1:num_DMs
                    C_temp(1,j+1) = sum(dot(U(:,j)'*resid_Y_temp, VY(j,:), 1));
                end
                C_temp(2:end,1) = C_temp(1,2:end)';
                for i=1:num_DMs
                    for j=1:num_DMs  
                        C_temp(i+1,j+1) = (U(:,i)'*U(:,j))*sum(dot(VY(i,:), VY(j,:), 1));
                    end
                end
                B_temp(1) = sum(dot(resid_Y_temp,X_temp,1));
                for k=1:num_DMs
                    B_temp(k+1) = sum(dot(U(:,k)*VY(k,:),X_temp,1));
                end

                % with autocorrelation term
                D = C_temp\B_temp;
                A_approx = D(1)*residual_matrix;
%                 for i = 1:num_DMs
%                     A_approx = A_approx + D(i+1) * U(:,i) * V(i,:);
%                 end
                A_approx = A_approx + U(:,1:num_DMs) * diag(D(2:end)) * V(1:num_DMs,:);
                A_approx = real(A_approx);
                loss = sum((X_temp_pred - A_approx * Y_temp_pred).^2,2);
                
                Y_test_plus_hat = A_approx  * Y_temp_pred;
                E_test = X_temp_pred - Y_test_plus_hat;                                     % Prediction error
                R2_DM = mean(1 - sum(E_test.^2, 2) ./ sum((X_temp_pred - mean(X_temp_pred, 2)).^2, 2));
                
                % null model - only consider autocorrelation
                c0 = sum(dot(Y_temp, Y_temp, 1));
                b0 = sum(dot(Y_temp,X_temp,1));
                a_null = real(b0/c0);
                loss_null = sum((X_temp_pred - a_null * Y_temp_pred).^2,2);
                Y_test_plus_hat = a_null  * Y_temp_pred;
                E_test = X_temp_pred - Y_test_plus_hat;                                     % Prediction error
                R2_null = mean(1 - sum(E_test.^2, 2) ./ sum((X_temp_pred - mean(X_temp_pred, 2)).^2, 2));
                
%                 loss_zero = sum((X_temp_pred - Y_temp_pred).^2,2);

                E_test = X_temp_pred - Y_temp_pred;                                     % Prediction error
                R2_zero = mean(1 - sum(E_test.^2, 2) ./ sum((X_temp_pred - mean(X_temp_pred, 2)).^2, 2));

                loss_session(n_ses) = mean(loss./loss_null);
                denom = sum((X_temp_pred - mean(X_temp_pred,2)).^2,2);
                R_squared_zero_session(n_ses) = R2_zero;
                R_squared_null_session(n_ses) = R2_null;
                R_squared_session(n_ses) = R2_DM;
            end
            
            loss_fold_current(n_test) = mean(loss_session);
            R_squared_fold_zero_current(n_test) = mean(R_squared_zero_session);
            R_squared_fold_null_current(n_test) = mean(R_squared_null_session);
            R_squared_fold_current(n_test) = mean(R_squared_session);
            session_count(n_test) = num_session;
        end
        loss_fold_current(isnan(loss_fold_current)) = 0;
        R_squared_fold_zero_current(isnan(R_squared_fold_zero_current)) = 0;
        R_squared_fold_null_current(isnan(R_squared_fold_null_current)) = 0;
        R_squared_fold_current(isnan(R_squared_fold_current)) = 0;
        loss_fold(n_cv,num_DMs/2) = sum(loss_fold_current.*session_count)/sum(session_count);
        R_squared_zero_fold(n_cv,num_DMs/2) = sum(R_squared_fold_zero_current.*session_count)/sum(session_count);
        R_squared_null_fold(n_cv,num_DMs/2) = sum(R_squared_fold_null_current.*session_count)/sum(session_count);
        R_squared_fold(n_cv,num_DMs/2) = sum(R_squared_fold_current.*session_count)/sum(session_count);
        disp(['The number of the included DM: ',num2str(num_DMs/2),', R_squared: ',num2str(R_squared_fold(n_cv,num_DMs/2),'%.6f'), ...
            ', R_squared (zero): ',num2str(R_squared_zero_fold(n_cv,num_DMs/2),'%.6f'),...
            ', R_squared (null): ',num2str(R_squared_null_fold(n_cv,num_DMs/2),'%.6f'),...
            ', Loss ratio compared to null model: ',num2str(loss_fold(n_cv,num_DMs/2),'%.6f')])
        toc
    end
    disp('===== Fold Done ====');
end

figure; 
subplot(2,1,1); bar(mean(R_squared_fold)-min(mean(R_squared_fold)));
subplot(2,1,2); bar(mean(loss_fold)-min(mean(loss_fold)));

save DMs/DM_cortical_subcortical_noROInorm_CV_prediction Phi_fold lambda_fold roi_exclude loss_fold R_squared_fold R_squared_zero_fold R_squared_null_fold

%% LASSO (training-test)
lambda_list =linspace(10,100,19);
%%% training %%%
R2_list_allsub_lambda = zeros(length(tau),length(lambda_list));
for nsub = 1:length(tau)
    disp(['=== (for training) sub-',num2str(nsub,'%04d'),' ===']);
    is_current_sub = data_index_sub == nsub;
    num_session = round(sum(is_current_sub)/(length(t_fine)-1));
    Y_sub = Y(:,is_current_sub);
    X_sub = X(:,is_current_sub);
    tic
    R2_list_ses_lambda = zeros(num_session,length(lambda_list));
    for n_ses = 1:num_session
        Y_sub_sess = Y_sub(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
        X_sub_sess = X_sub(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
        
        n_inner_cv = 5;
        Y_fitting = Y_sub_sess(:,1:round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv));
        X_fitting = X_sub_sess(:,1:round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv));
        Y_validation = Y_sub_sess(:,round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv)+1:end);
        X_validation = X_sub_sess(:,round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv)+1:end);
        
        R2_list_lambda = linear_LASSO_yj_training(Y_fitting,X_fitting,Y_validation,X_validation,lambda_list);
        [best_R2,best_idx] = max(R2_list_lambda);
        fprintf('Current best h = %.3f, best R2 = %.6f \n', lambda_list(best_idx), best_R2);
        R2_list_ses_lambda(n_ses,:) = R2_list_lambda;
    end
    toc
    R2_list_allsub_lambda(nsub,:) = mean(R2_list_ses_lambda);
end

save results/LASSO_training_temp R2_list_allsub_lambda lambda_list

%%% test %%%
R_squared_fold_LASSO = zeros(cv_num,1);
for n_cv = 1:cv_num
    disp(['=== CV prediction fold #',num2str(n_cv),'===']);
    lambda = lambda_fold{n_cv};
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    temp_R2_lambda_list = R2_list_allsub_lambda(current_train_subjects,:);
    temp_tau = tau(current_train_subjects);
    temp_R2_lambda_list(temp_tau==0,:) = [];
    temp_lambda_list = lambda_list;
    
    temp_lambda_list(1) = [];
    temp_R2_lambda_list(:,1) = [];
    
    mean_temp_R2_lambda_list = mean(temp_R2_lambda_list);
    lambda_fine = temp_lambda_list(1):0.001:temp_lambda_list(end);
    mean_temp_R2_lambda_fine = spline(temp_lambda_list,mean_temp_R2_lambda_list,lambda_fine);
    
    [~,idx_max] = max(mean_temp_R2_lambda_fine);
    lambda_max_training = lambda_fine(idx_max);
    fprintf('Hyperparameter lambda of current fold: %.3f\n',lambda_max_training);
    
    R_squared_fold_current = zeros(1,length(current_test_subjects));
    session_count = zeros(1,length(current_test_subjects));
    for n_test = 1:length(current_test_subjects)
        tic
        current_subject = current_test_subjects(n_test);
        is_test = data_index_sub == current_subject;
        num_session = round(sum(is_test)/(length(t_fine)-1));
        X_test = X(:,is_test);
        Y_test = Y(:,is_test);
        R_squared_session = zeros(1,num_session);
        for n_ses = 1:num_session
            X_test_sess = X_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
            Y_test_sess = Y_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));

            X_temp = X_test_sess(:,1:round(length(t_fine)*4/5));
            Y_temp = Y_test_sess(:,1:round(length(t_fine)*4/5));
            X_temp_pred = X_test_sess(:,round(length(t_fine)*4/5)+1:end);
            Y_temp_pred = Y_test_sess(:,round(length(t_fine)*4/5)+1:end);
            
            R2 = linear_LASSO_yj(Y_temp,X_temp,Y_temp_pred,X_temp_pred,lambda_max_training);
            R_squared_session(n_ses) = R2;
        end
        R_squared_fold_current(n_test) = mean(R_squared_session);
        session_count(n_test) = num_session;
        fprintf('R-squared (current test subject):%.6f \n', mean(R_squared_session));
        toc
        disp('');
    end
    R_squared_fold_current(isnan(R_squared_fold_current)) = 0;
    R_squared_fold_LASSO(n_cv) = sum(R_squared_fold_current.*session_count)/sum(session_count);
    disp(['R_squared (current fold): ',num2str(R_squared_fold_LASSO(n_cv),'%.6f'),'\n\n'])
end

save('results/other_model_cv_results', 'R_squared_fold_LASSO', '-append');

%% Hyperparameter Plot for All CV Folds in a Single Figure

% Determine the number of CV folds
% Ensure cv_num is defined; if not, set it based on the length of lambda_fold
if ~exist('cv_num', 'var') || isempty(cv_num)
    cv_num = length(lambda_fold);
end

% Determine the layout for subplots (rows and columns)
% For better visualization, aim for a nearly square grid
num_cols = ceil(sqrt(cv_num));
num_rows = ceil(cv_num / num_cols);

% Create a new figure for all CV folds
figure('Name', 'R² vs \lambda Across CV Folds', 'NumberTitle', 'off','Position',[100, 100, 1400, 700]);

% Create a tiled layout
t = tiledlayout(num_rows, num_cols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Initialize a cell to store optimal lambdas for all folds (optional)
optimal_lambdas = zeros(cv_num, 1);

% Loop through each CV fold
for n_cv = 1:cv_num
    disp(['=== CV Prediction Fold #', num2str(n_cv), ' ===']);
    lambda = lambda_fold{n_cv};
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    temp_R2_lambda_list = R2_list_allsub_lambda(current_train_subjects, :);
    temp_tau = tau(current_train_subjects);
    temp_R2_lambda_list(temp_tau == 0, :) = [];
    temp_lambda_list = lambda_list;
    
    % Remove the first lambda value and corresponding R2 values
    temp_lambda_list(1) = [];
    temp_R2_lambda_list(:,1) = [];
    
    % Calculate the mean R2 across subjects for each lambda
    mean_temp_R2_lambda_list = mean(temp_R2_lambda_list);
    
    % Create a finer lambda grid for interpolation
    lambda_fine = temp_lambda_list(1):0.001:temp_lambda_list(end);
    
    % Perform spline interpolation on the mean R2 values
    mean_temp_R2_lambda_fine = spline(temp_lambda_list, mean_temp_R2_lambda_list, lambda_fine);
    
    % Find the lambda that maximizes the interpolated mean R2
    [~, idx_max] = max(mean_temp_R2_lambda_fine);
    lambda_max_training = lambda_fine(idx_max);
    max_R2 = mean_temp_R2_lambda_fine(idx_max);
    
    fprintf('Optimal lambda for CV fold %d: %.3f\n', n_cv, lambda_max_training);
    
    % Store the optimal lambda (optional)
    optimal_lambdas(n_cv) = lambda_max_training;
    
    %%% Plotting within the tiled layout
    % Select the next tile
    nexttile;
    hold on; % Hold on to plot multiple elements
    
    % Plot the original mean R2 values as blue dots
    plot(temp_lambda_list, mean_temp_R2_lambda_list, 'bo', ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Original R²');
    
    % Define sky blue color using RGB triplet
    skyblue = [0.529, 0.808, 0.980];
    % Plot the interpolated mean R2 values as a sky blue line
    plot(lambda_fine, mean_temp_R2_lambda_fine, '-', ...
        'Color', skyblue, 'LineWidth', 1.5, 'DisplayName', 'Interpolated R²');
    
    % Plot the optimal lambda as a red 'x'
    plot(lambda_max_training, max_R2, 'rx', ...
        'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Optimal \lambda');
    
     % -------------------------------
    % Add Text Annotation for Optimal Lambda
    % -------------------------------
    % Define the text label string with lambda value formatted to three decimal places
    lambda_text = sprintf('\\lambda_{best} = %.3f', lambda_max_training);
    
    % Determine text position offsets for better visibility
    % Adjust the offsets as needed based on your data range
    x_offset = 0.05 * (max(temp_lambda_list) - min(temp_lambda_list));
    y_offset = 0.3 * (max(mean_temp_R2_lambda_fine) - min(mean_temp_R2_lambda_fine));
    
    % Add text near the optimal lambda point
    text(lambda_max_training + x_offset, max_R2 - y_offset, lambda_text, ...
        'Color', 'black', 'FontSize', 18, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    
    hold off; % Release the hold
    
    % Enhance subplot aesthetics
    xlabel('\lambda', 'FontSize', 18);
    ylabel('Mean R²', 'FontSize', 18);
    title(['Fold ', num2str(n_cv)], 'FontSize', 22);
    if n_cv == cv_num
        legend('Location', 'best', 'FontSize', 18);
    end
    grid on; % Add grid lines for better readability
end

% Adjust the overall figure title
title(t, 'R² vs \lambda Across All CV Folds', 'FontSize', 16, 'FontWeight', 'bold');

% Optionally, display all optimal lambdas after the plots
fprintf('\nOptimal \lambda values across all CV folds:\n');
disp(optimal_lambdas);



%% Non-linear (manifold)
% h_list = [10.02,11.54];
h_list = 2:2:20;
%%% training %%%
R2_list_allsub_h = zeros(length(tau),length(h_list));
for nsub = 1:length(tau)
    disp(['=== (for training) sub-',num2str(nsub,'%04d'),' ===']);
    is_current_sub = data_index_sub == nsub;
    num_session = round(sum(is_current_sub)/(length(t_fine)-1));
    Y_sub = Y(:,is_current_sub);
    X_sub = X(:,is_current_sub);
    tic
    R2_list_ses_h = zeros(num_session,length(h_list));
    for n_ses = 1:num_session
        Y_sub_sess = Y_sub(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
        X_sub_sess = X_sub(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
        
        n_inner_cv = 5;
        Y_fitting = Y_sub_sess(:,1:round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv));
        X_fitting = X_sub_sess(:,1:round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv));
        Y_validation = Y_sub_sess(:,round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv)+1:end);
        X_validation = X_sub_sess(:,round(size(Y_sub_sess,2)*(n_inner_cv-1)/n_inner_cv)+1:end);
        
        R2_list_h = zeros(1,length(h_list));
        for n_h = 1:length(h_list)
            h = h_list(n_h);
            R2 = nonlinear_manifold_yj(Y_fitting,X_fitting,Y_validation,X_validation,h);
            R2_list_h(n_h) = R2;
        end
        [best_R2,best_idx] = max(R2_list_h);
        fprintf('Current best h = %.3f, best R2 = %.6f \n', h_list(best_idx), best_R2);
        R2_list_ses_h(n_ses,:) = R2_list_h;
    end
    toc
    R2_list_allsub_h(nsub,:) = mean(R2_list_ses_h);
end

save results/manifold_training_temp R2_list_allsub_h h_list

%%% test %%%
R_squared_fold_manifold = zeros(cv_num,1);
for n_cv = 1:cv_num
    disp(['=== CV prediction fold #',num2str(n_cv),'===']);
    lambda = lambda_fold{n_cv};
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    temp_R2_h_list = R2_list_allsub_h(current_train_subjects,:);
    temp_tau = tau(current_train_subjects);
    temp_R2_h_list(temp_tau==0,:) = [];
    temp_h_list = h_list;
    
%     temp_h_list(1) = [];
%     temp_R2_h_list(:,1) = [];
    
    mean_temp_R2_h_list = mean(temp_R2_h_list);
    h_fine = temp_h_list(1):0.001:temp_h_list(end);
    mean_temp_R2_h_fine = spline(temp_h_list,mean_temp_R2_h_list,h_fine);
    
    [~,idx_max] = max(mean_temp_R2_h_fine);
    h_max_training = h_fine(idx_max);
    fprintf('Hyperparameter h of current fold: %.3f\n',h_max_training);
    
    R_squared_fold_current = zeros(1,length(current_test_subjects));
    session_count = zeros(1,length(current_test_subjects));
    for n_test = 1:length(current_test_subjects)
        tic
        current_subject = current_test_subjects(n_test);
        is_test = data_index_sub == current_subject;
        num_session = round(sum(is_test)/(length(t_fine)-1));
        X_test = X(:,is_test);
        Y_test = Y(:,is_test);
        R_squared_session = zeros(1,num_session);
        for n_ses = 1:num_session
            X_test_sess = X_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
            Y_test_sess = Y_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));

            X_temp = X_test_sess(:,1:round(length(t_fine)*4/5));
            Y_temp = Y_test_sess(:,1:round(length(t_fine)*4/5));
            X_temp_pred = X_test_sess(:,round(length(t_fine)*4/5)+1:end);
            Y_temp_pred = Y_test_sess(:,round(length(t_fine)*4/5)+1:end);
            
            R2 = nonlinear_manifold_yj(Y_temp,X_temp,Y_temp_pred,X_temp_pred,h_max_training);
            R_squared_session(n_ses) = R2;
        end
        R_squared_fold_current(n_test) = mean(R_squared_session);
        session_count(n_test) = num_session;
        fprintf('R-squared (current test subject):%.6f \n', mean(R_squared_session));
        toc
        disp('');
    end
    R_squared_fold_current(isnan(R_squared_fold_current)) = 0;
    R_squared_fold_manifold(n_cv) = sum(R_squared_fold_current.*session_count)/sum(session_count);
    disp(['R_squared (current fold): ',num2str(R_squared_fold_manifold(n_cv),'%.6f'),'\n\n'])
end

save('results/other_model_cv_results', 'R_squared_fold_manifold', '-append');

%% Hyperparameter Plot for All CV Folds in a Single Figure (Using Hyperparameter h)

% Determine the number of CV folds
% Ensure cv_num is defined; if not, set it based on the length of h_fold
if ~exist('cv_num', 'var') || isempty(cv_num)
    cv_num = length(h_fold);
end

% Determine the layout for subplots (rows and columns)
% For better visualization, aim for a nearly square grid
num_cols = ceil(sqrt(cv_num));
num_rows = ceil(cv_num / num_cols);

% Create a new figure for all CV folds
figure('Name', 'R² vs h Across CV Folds', 'NumberTitle', 'off', 'Position', [100, 100, 1400, 700]);

% Create a tiled layout
t = tiledlayout(num_rows, num_cols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Initialize a vector to store optimal h values for all folds (optional)
optimal_h = zeros(cv_num, 1);

% Loop through each CV fold
for n_cv = 1:cv_num
    disp(['=== CV Prediction Fold #', num2str(n_cv), ' ===']);
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    temp_R2_h_list = R2_list_allsub_h(current_train_subjects, :);
    temp_tau = tau(current_train_subjects);
    temp_R2_h_list(temp_tau == 0, :) = [];
    temp_h_list = h_list;
    
    % Remove the first h value and corresponding R² values
    temp_h_list(1) = [];
    temp_R2_h_list(:,1) = [];
    
    % Calculate the mean R² across subjects for each h
    mean_temp_R2_h_list = mean(temp_R2_h_list);
    
    % Create a finer h grid for interpolation
    h_fine = temp_h_list(1):0.001:h_list(end);
    
    % Perform spline interpolation on the mean R² values
    mean_temp_R2_h_fine = spline(temp_h_list, mean_temp_R2_h_list, h_fine);
    
    % Find the h that maximizes the interpolated mean R²
    [~, idx_max] = max(mean_temp_R2_h_fine);
    h_max_training = h_fine(idx_max);
    max_R2 = mean_temp_R2_h_fine(idx_max);
    
    fprintf('Optimal h for CV fold %d: %.3f\n', n_cv, h_max_training);
    
    % Store the optimal h (optional)
    optimal_h(n_cv) = h_max_training;
    
    %%% Plotting within the tiled layout
    % Select the next tile
    nexttile;
    hold on; % Hold on to plot multiple elements
    
    % Plot the original mean R² values as blue dots
    plot(temp_h_list, mean_temp_R2_h_list, 'bo', ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Original R²');
    
    % Define sky blue color using RGB triplet
    skyblue = [0.529, 0.808, 0.980];
    % Plot the interpolated mean R² values as a sky blue line
    plot(h_fine, mean_temp_R2_h_fine, '-', ...
        'Color', skyblue, 'LineWidth', 1.5, 'DisplayName', 'Interpolated R²');
    
    % Plot the optimal h as a red 'x'
    plot(h_max_training, max_R2, 'rx', ...
        'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Optimal h');
    
    % -------------------------------
    % Add Text Annotation for Optimal h
    % -------------------------------
    % Define the text label string with h value formatted to three decimal places
    h_text = sprintf('h_{best} = %.3f', h_max_training);
    
    % Determine text position offsets for better visibility
    % Adjust the offsets as needed based on your data range
    x_offset = 0.05 * (max(temp_h_list) - min(temp_h_list));
    y_offset = 0.3 * (max(mean_temp_R2_h_fine) - min(mean_temp_R2_h_fine));
    
    % Add text near the optimal h point
    text(h_max_training + x_offset, max_R2 - y_offset, h_text, ...
        'Color', 'black', 'FontSize', 18, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    
    hold off; % Release the hold
    
    % Enhance subplot aesthetics
    xlabel('h', 'FontSize', 18);
    ylabel('Mean R²', 'FontSize', 18);
    title(['Fold ', num2str(n_cv)], 'FontSize', 22);
    if n_cv == cv_num
        legend('Location', 'best', 'FontSize', 18);
    end
    grid on; % Add grid lines for better readability
end

% Adjust the overall figure title
title(t, 'R² vs h Across All CV Folds', 'FontSize', 16, 'FontWeight', 'bold');

% Optionally, display all optimal h values after the plots
fprintf('\nOptimal h values across all CV folds:\n');
disp(optimal_h);

