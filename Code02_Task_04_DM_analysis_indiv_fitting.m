clear; clc;
current_path = pwd;

task_names = {'WM','EMOTION','MOTOR','LANGUAGE','GAMBLING','SOCIAL','RELATIONAL'};

for n_task = 1:length(task_names)
    task = task_names{n_task};

    load(['results/HCP_timeseries_tfMRI_',task,'_cortical_subcortical_extracted_filtered_CompCor_meta.mat']);
    load(['results/HCP_timeseries_tfMRI_',task,'_cortical_subcortical_extracted_filtered_CompCor.mat']);

    %%
    n_time = len_time;

    t_sample = 0.72;
    TRtarget = 0.72;
%     TRtarget = 1.5;

    t = (1:n_time) * (t_sample);
    t_fine = TRtarget:TRtarget:t(end);

    num_sub = size(time_series_preproc_filtered,1);
    tau=zeros(num_sub,1);

    i_num = 0;
    sub_ids_estimated = cell(0);
    %%
    load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm

    num_DMs = 10;

    N = 716;

    U=Phi_sorted(:,1:num_DMs);
    V=pinv(U);
    residual_matrix = eye(size(U,1)) - U*V;

    D = zeros(num_DMs+1,size(time_series_preproc_filtered,1));
    loss_DM_model = zeros(1,size(time_series_preproc_filtered,1));
    for nsub = 1:size(time_series_preproc_filtered,1)
        tau_temp = 0;
        disp(['start: sub#' num2str(nsub)]);
        X_temp = []; Y_temp = []; C_temp = [];
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
            C_temp = zeros(num_DMs+1,num_DMs+1);
            B_temp = zeros(num_DMs+1,1);

            resid_Y_temp = residual_matrix * Y_temp;
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

            D(:,nsub) = C_temp\B_temp;
            disp(['end: sub#' num2str(nsub)]);


        end
        cd(current_path);
        tau(nsub) = tau_temp;
        toc
    end


    save(['DMs/DM_tfMRI_',task,'_cortical_subcortical_ext_fbDMD_noROInorm_indiv_10'],'tau', 'D', 'sub_ids')

end