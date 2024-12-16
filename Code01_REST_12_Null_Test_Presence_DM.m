clear; clc;

load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm
%% Generate random data with band-pass filtering. Subject-level DM fitting
num_Perm = 100;

num_DMs = 10;

U=Phi_sorted(:,1:num_DMs);
V=pinv(U);
residual_matrix = eye(size(U,1)) - U*V;
D=zeros(num_DMs+1,num_Perm);

for nperm = 1:num_Perm
    disp(nperm);
    
    random_data = randn(716,1500);
    %%% band pass filter
    % Define the sampling frequency
    Fs = 1/1.5;  % Hz

    % Define the frequency range for the band-pass filter
    lowCutoff = 0.01;  % Hz
    highCutoff = 0.1;  % Hz
    
    % Apply the filter
    filteredData = bandpass(random_data',[lowCutoff, highCutoff],Fs)';
    
    
    Z_temp = zeros(num_DMs+1,num_DMs+1);
    W_temp = zeros(num_DMs+1,1);
    
    X_temp = filteredData(:,2:end);
    Y_temp = filteredData(:,1:end-1);
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

    D(:,nperm) = Z_temp\W_temp;
end

%% Statistical test
disp('##### Statistical Test #####')
disp('Since no DM present, t-test using angles should not be significant for all DMs ...');
angel_D = angle(D(2:end,:));
for n_dm = 1:2:num_DMs
    [h,p] = ttest(angel_D(n_dm,:));
    disp(['DM #',num2str((n_dm+1)/2), ': p = ', num2str(p)])
end