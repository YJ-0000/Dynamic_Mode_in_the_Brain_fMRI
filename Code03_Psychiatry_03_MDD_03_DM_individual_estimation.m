clear; clc;
current_path = pwd;

load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm
load results/MDD_timeseires_cortical_subcortical_denoised_filtered
%%
count_roi = zeros(size(Phi_sorted,1),1);
for ndir = 1:length(MDD_timeseries_extracted_denoised_filtered)
    temp_timeseries_extracted_denoised_filtered = MDD_timeseries_extracted_denoised_filtered{ndir};
    for nsub = 1:length(temp_timeseries_extracted_denoised_filtered)
        y = temp_timeseries_extracted_denoised_filtered{nsub};
        y(roi_exclude,:) = [];
        var_y = var(y,0,2);
        count_roi = count_roi + (var_y<0.0001);
    end
end

%%
TRtarget = 1.5;
num_DMs = 10;
thres = 1;

U=Phi_sorted(:,1:num_DMs);
U(count_roi>thres,:)=[];
V=pinv(U);
residual_matrix = eye(size(U,1)) - U*V;

D_list = cell(size(MDD_timeseries_extracted_denoised_filtered));
B_list = cell(size(MDD_timeseries_extracted_denoised_filtered));
for ndir = 1:length(MDD_timeseries_extracted_denoised_filtered) 
    temp_timeseries_extracted_denoised_filtered = MDD_timeseries_extracted_denoised_filtered{ndir};
    temp_TR_info = MDD_TR_info{ndir};
    temp_sub_list = sub_list{ndir};
    
    D_temp = zeros(num_DMs+1,length(temp_timeseries_extracted_denoised_filtered));
    B_temp = zeros(num_DMs,length(temp_timeseries_extracted_denoised_filtered));
    for nsub = 1:length(temp_timeseries_extracted_denoised_filtered)
        tr = temp_TR_info(nsub);
        sub_name = temp_sub_list{nsub};
        fprintf(['sub >> ',sub_name,'.\n']);
        
        
        y = temp_timeseries_extracted_denoised_filtered{nsub};
        Z_temp = zeros(num_DMs+1,num_DMs+1);
        W_temp = zeros(num_DMs+1,1);

        tic
        
        y(roi_exclude,:) = [];
        y(count_roi>thres,:)=[];
        t = (1:size(y,2)) * (tr);
        t_fine = TRtarget:TRtarget:t(end);
        if tr ~= TRtarget
            pp = spline(t, y);          % Compute piecewise polynomial (B-spline) representation
            y_fine = ppval(pp, t_fine); % Evaluate the piecewise polynomial at the finer time points
        else
            y_fine = y;
        end

        X_temp = y_fine(:,2:end);
        Y_temp = y_fine(:,1:end-1);
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
        toc

        D_temp(:,nsub) = Z_temp\W_temp;
        B_temp(:,nsub) = mean(abs(VY),2);
        
    end
    D_list{ndir} = D_temp;
    B_list{ndir} = B_temp;
end

%%
save DMs/DM_MDD_cortical_subcortical_noROInorm_indiv_10 D_list B_list Phi_sorted lambda count_roi thres sub_list sub_info