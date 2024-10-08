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
    lambda(idx_exclude) = [];
    Phi_sorted(:,idx_exclude) = [];
    [lambda,idx_sort] = sort(lambda,'descend');
    Phi_sorted = Phi_sorted(:,idx_sort);
    
    toc
    disp(['=== DONE ===']);
    
    lambda_fold{n_cv} = lambda;
    Phi_fold{n_cv} = Phi_sorted;
end

%% save CV results
save DMs/DM_cortical_subcortical_noROInorm_CV Phi_fold lambda_fold roi_exclude

%% CV prediction
disp('*** CV prediction started! ***');
max_DMs = 24;
% loss_fold = zeros(cv_num,max_DMs/2);
for n_cv = 2:cv_num
    disp(['=== CV prediction fold #',num2str(n_cv),'===']);
    lambda = lambda_fold{n_cv};
    Phi_sorted = Phi_fold{n_cv};
    
    current_train_subjects = train_id_list{n_cv};
    current_test_subjects = test_id_list{n_cv};
    
    for num_DMs = 2:2:max_DMs
        tic
        U=Phi_sorted(:,1:num_DMs);
        V=pinv(U);
        
        loss_fold_current = zeros(1,length(current_test_subjects));
        session_count = zeros(1,length(current_test_subjects));
        for n_test = 1:length(current_test_subjects)
            current_subject = current_test_subjects(n_test);
            is_test = data_index_sub == current_subject;
%             disp(['Data points in X_train: ', num2str(sum(is_test))]);

            num_session = round(sum(is_test)/(length(t_fine)-1));
            X_test = X(:,is_test);
            Y_test = Y(:,is_test);
            loss_session = zeros(1,num_session);
            for n_ses = 1:num_session
                X_test_sess = X_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
                Y_test_sess = Y_test(:,(n_ses-1)*(length(t_fine)-1)+1:n_ses*(length(t_fine)-1));
                C_temp = zeros(num_DMs+1,num_DMs+1);
                B_temp = zeros(num_DMs+1,1);

                X_temp = X_test_sess(:,1:round(length(t_fine)*4/5));
                Y_temp = Y_test_sess(:,1:round(length(t_fine)*4/5));
                X_temp_pred = X_test_sess(:,round(length(t_fine)*4/5)+1:end);
                Y_temp_pred = Y_test_sess(:,round(length(t_fine)*4/5)+1:end);
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


                D = C_temp\B_temp;

                A_approx = D(1)*eye(size(X_temp,1));
                for i = 1:num_DMs
                    A_approx = A_approx + D(i+1) * U(:,i) * V(i,:);
                end
                A_approx = real(A_approx);
                loss = sum((X_temp_pred - A_approx * Y_temp_pred).^2,'all');

                a_null = real(B_temp(1)/C_temp(1));
                loss_null = sum((X_temp_pred - a_null * Y_temp_pred).^2,'all');

                loss_session(n_ses) = loss/loss_null;
            end
            
            loss_fold_current(n_test) = mean(loss_session);
            session_count(n_test) = num_session;
        end
        loss_fold_current(isnan(loss_fold_current)) = 0;
        loss_fold(n_cv,num_DMs/2) = sum(loss_fold_current.*session_count)/sum(session_count);
        toc
    end
    disp('===== Fold Done ====');
end

figure; bar(mean(loss_fold)-min(mean(loss_fold)));

save DMs/DM_cortical_subcortical_noROInorm_CV_prediction Phi_fold lambda_fold roi_exclude loss_fold

