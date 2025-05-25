clear; clc;
current_path = pwd;

task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};

for n_task = 1:length(task_names)
    task = task_names{n_task};

    load(['results/HCP_timeseries_tfMRI_',task,'_cortical_subcortical_extracted_filtered_CompCor_meta.mat']);
    load(['results/HCP_timeseries_tfMRI_',task,'_cortical_subcortical_extracted_filtered_CompCor.mat']);
    load(['results/HCP_timeseries_tfMRI_',task,'_cortical_subcortical_extracted_task_input']);
    load(['results/HCP_timeseries_tfMRI_',task,'_subject_exclude_info']);
    
    %%
    is_sub_exclude = true;
    if is_sub_exclude
        for nsub = 1:length(sub_ids)
            if does_have_MMSE(nsub) || is_cognitive_impaired(nsub) || is_RL_processing_errors(nsub)
                time_series_preproc_filtered(nsub,:) = {[],[]}; %#ok<SAGROW>
            else
                if is_excluded_due_movement(nsub,1)
                    time_series_preproc_filtered(nsub,1:2) = {[],[]}; %#ok<SAGROW>
                end
            end
        end
    end

    remaining_sub_idx = false(length(sub_ids),1);
    for nsub = 1:length(sub_ids)
        if ~isempty(time_series_preproc_filtered{nsub,1}) || ~isempty(time_series_preproc_filtered{nsub,2}) 
            remaining_sub_idx(nsub) = true;
        end
    end

    load secure_data/path_info;
    gene_data_table = readtable(gene_data_path,'VariableNamingRule','preserve');
    behav_data_table = readtable(behav_data_path,'VariableNamingRule','preserve');
    for nrow = size(gene_data_table,1):-1:1
        if ~any(sub_ids==gene_data_table(nrow,'Subject').Variables)
            gene_data_table(nrow,:) = [];
        end
    end
    for nrow = size(behav_data_table,1):-1:1
        if ~any(sub_ids==behav_data_table(nrow,'Subject').Variables)
            behav_data_table(nrow,:) = [];
        end
    end
    gene_data_table = sortrows(gene_data_table, 'Subject');
    behav_data_table = sortrows(behav_data_table, 'Subject');

    ages = gene_data_table.Age_in_Yrs;
    genders = behav_data_table.Gender;

    num_female = sum(strcmp(genders(remaining_sub_idx),'F'));
    mean_age = mean(ages(remaining_sub_idx));
    std_age = std(ages(remaining_sub_idx));

    fprintf('Total number of subjects: %d, Female=%d, mean age=%0.2f, std=%0.2f \n', ...
        sum(remaining_sub_idx),num_female,mean_age,std_age);

    %%
    n_time = len_time;

    t_sample = 0.72;
%     TRtarget = 0.72;
    TRtarget = 1.5;

    t = (1:n_time) * (t_sample);
    t_fine = TRtarget:TRtarget:t(end);

    num_sub = size(time_series_preproc_filtered,1);
    tau=zeros(num_sub,1);

    i_num = 0;
    sub_ids_estimated = cell(0);
    %%
    load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_subExclude

    num_DMs = 10;

    N = 716;

    U=Phi_sorted(:,1:num_DMs);
    V=pinv(U);
    residual_matrix = eye(size(U,1)) - U*V;

    D = zeros(num_DMs+1,size(time_series_preproc_filtered,1));
    B = zeros(num_DMs,size(time_series_preproc_filtered,1));
    loss_DM_model = zeros(1,size(time_series_preproc_filtered,1));
    for nsub = 1:size(time_series_preproc_filtered,1)
        tau_temp = 0;
        disp(['start: sub#' num2str(nsub)]);
        X_temp = []; Y_temp = []; Z_temp = [];
        tic
        for nses = 1:2
            if ~isempty(time_series_preproc_filtered{nsub,nses})
                i_num = i_num + 1;
                y = time_series_preproc_filtered{nsub,nses};
                if t_sample ~= TRtarget
                    pp = spline(t, y);          % Compute piecewise polynomial (B-spline) representation
                    y_fine = ppval(pp, t_fine); % Evaluate the piecewise polynomial at the finer time points
                else
                    y_fine = y;
                end
                tau_temp = tau_temp + (size(y_fine,2)-1);

                y_fine(roi_exclude,:) = [];

                X_temp = [X_temp,y_fine(:,2:end)];
                Y_temp = [Y_temp,y_fine(:,1:end-1)];
            end
        end

        if tau_temp > 0
            Z_temp = zeros(num_DMs+1,num_DMs+1);
            W_temp = zeros(num_DMs+1,1);

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

            D(:,nsub) = Z_temp\W_temp;
            B(:,nsub) = mean(abs(VY),2);
            disp(['end: sub#' num2str(nsub)]);


        end
        cd(current_path);
        tau(nsub) = tau_temp;
        toc
    end


    save(['DMs/DM_tfMRI_',task,'_cortical_subcortical_ext_fbDMD_noROInorm_subExclude_indiv_10_B'],'tau', 'D', 'B', 'sub_ids','TRtarget','remaining_sub_idx')

end