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
    idx_exclude = (abs(angle(lambda)) < 2*pi*1.5*0.01) | (abs(lambda)>1);
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
max_DMs = 18;

Phi_abs_all = [];
for n_cv = 1:cv_num
    Phi = Phi_fold{n_cv};
    Phi_abs_all = [Phi_abs_all,abs(Phi(:,1:2:max_DMs))]; %#ok<AGROW>
end
[cluster_idx, C] = kmeans(Phi_abs_all', max_DMs/2, 'Replicates', 10);
for n_cv = 1:cv_num
    assert(length(unique(cluster_idx((n_cv-1)*max_DMs/2+1:n_cv*max_DMs/2))) == max_DMs/2)
end

consistency_list = cell(max_DMs/2,1);
for n_dm = 1:max_DMs/2
    figure;
    iter_n = 0;
    Phi_select = Phi_abs_all(:,cluster_idx==n_dm);
    corr_mat = zeros(cv_num);
    for n_cv1 = 1:cv_num
        for n_cv2 = 1:cv_num
            iter_n = iter_n + 1;
            subplot(cv_num,cv_num,iter_n);
            scatter(Phi_select(:,n_cv1),Phi_select(:,n_cv2))
            
            r = corrcoef(Phi_select(:,n_cv1),Phi_select(:,n_cv2));
            corr_mat(n_cv1,n_cv2) = r(2);
        end
    end
    consistency_list{n_dm} = corr_mat;
end

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
                
                % null model - only consider autocorrelation
                c0 = sum(dot(Y_temp, Y_temp, 1));
                b0 = sum(dot(Y_temp,X_temp,1));
                a_null = real(b0/c0);
                loss_null = sum((X_temp_pred - a_null * Y_temp_pred).^2,2);
                
                loss_zero = sum((X_temp_pred - Y_temp_pred).^2,2);

                loss_session(n_ses) = mean(loss./loss_null);
                denom = sum((X_temp_pred - mean(X_temp_pred,2)).^2,2);
                R_squared_zero_session(n_ses) = 1-mean(loss_zero./denom);
                R_squared_null_session(n_ses) = 1-mean(loss_null./denom);
                R_squared_session(n_ses) = 1-mean(loss./denom);
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

%% For comparison purpose

R_squared_fold_lasso = zeros(cv_num,1);
for n_cv = 1:cv_num
    disp(['=== CV prediction fold #',num2str(n_cv),'===']);
    lambda = lambda_fold{n_cv};
    Phi_sorted = Phi_fold{n_cv};
    Phi_all = Phi_all_fold{n_cv};
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    R_squared_fold_current = zeros(1,length(current_test_subjects));
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
            
            Y_temp_diff = X_temp - Y_temp;
            
            % Initialize Theta for multiple response variables
            lambda_list =linspace(10,100,19);
            [num_responses, ~] = size(Y_temp_diff);
            Theta = zeros(num_responses, size(Y_temp, 1),length(lambda_list));
            
            MSE_accum = zeros(length(lambda_list),1);
            
            n_inner_cv = 4;
            % Perform Lasso regression for each response variable
            parfor i = 1:num_responses
                % Transpose data to match lasso's expected input
                [B, ~] = lasso(Y_temp(:,1:round(size(Y_temp,2)*(n_inner_cv-1)/n_inner_cv))', ...
                    Y_temp_diff(i, 1:round(size(Y_temp,2)*(n_inner_cv-1)/n_inner_cv))', 'lambda', lambda_list, 'Standardize', true);
                % If multiple coefficients are returned, select the one corresponding to lambda
                Theta(i, :, :) = B;  % Assuming lambda is a scalar
            end
            for n_inner = 1:length(lambda_list)
                W = (eye(size(Theta,1)) + squeeze(Theta(:,:,n_inner)));
                Y_hat = W  * Y_temp(:,round(size(Y_temp,2)*(n_inner_cv-1)/n_inner_cv)+1:end);
                loss = sum((X_temp(:,round(size(Y_temp,2)*(n_inner_cv-1)/n_inner_cv)+1:end) - Y_hat).^2,2);
                denom = sum((X_temp(:,round(size(Y_temp,2)*(n_inner_cv-1)/n_inner_cv)+1:end) - mean(X_temp(:,round(size(Y_temp,2)*(n_inner_cv-1)/n_inner_cv)+1:end),2)).^2,2);
                R_squared_temp = 1-mean(loss./denom);
                MSE_accum(n_inner) = R_squared_temp;
            end
            [~,max_idx] = max(MSE_accum);
            best_lambda = lambda_list(max_idx);
            fprintf('Current best hyperpameter = %.f \n',best_lambda);
            
            Theta_best = zeros(num_responses, size(Y_temp, 1));
            parfor i = 1:num_responses
                % Transpose data to match lasso's expected input
                [B, ~] = lasso(Y_temp', Y_temp_diff(i, :)', 'lambda', best_lambda, 'Standardize', true);
                % If multiple coefficients are returned, select the one corresponding to lambda
                Theta_best(i, :) = B;  % Assuming lambda is a scalar
            end
            
            W = (eye(size(Theta_best)) + Theta_best) ;
            Y_hat = W  * Y_temp_pred;
            loss = sum((X_temp_pred - Y_hat).^2,2);
            denom = sum((X_temp_pred - mean(X_temp_pred,2)).^2,2);
            R_squared_session(n_ses) = 1-mean(loss./denom);
        end
        R_squared_fold_current(n_test) = mean(R_squared_session);
        disp(mean(R_squared_session));
        toc
    end
    R_squared_fold_current(isnan(R_squared_fold_current)) = 0;
    R_squared_fold_lasso(n_cv) = sum(R_squared_fold_current.*session_count)/sum(session_count);
    disp(['R_squared: ',num2str(R_squared_fold_lasso(n_cv),'%.6f')])
end

save results/other_model_cv_results R_squared_fold_lasso