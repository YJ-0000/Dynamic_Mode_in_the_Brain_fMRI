clear; clc;
current_path = pwd;
conn;
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

            tau_temp = tau_temp + (size(y_fine,2)-1);
            X(:,(i_num-1)*(length(t_fine)-1)+1:i_num*(length(t_fine)-1)) = y_fine(:,2:end);
            Y(:,(i_num-1)*(length(t_fine)-1)+1:i_num*(length(t_fine)-1)) = y_fine(:,1:end-1);
        end
    end
    tau(nsub) = tau_temp;
end
X(:,i_num*(length(t_fine)-1)+1:end) = [];
Y(:,i_num*(length(t_fine)-1)+1:end) = [];

disp('*** total number of time points ***');
disp([i_num*(length(t_fine)-1), size(X,2)]);

% removing rois with no signal
var_X = var(X,0,2);
roi_exclude = var_X < 0.001;
X(roi_exclude,:) = [];
Y(roi_exclude,:) = [];

clear time_series_denoised_filtered

%%
%%% fbDMD
disp('*** Extended fbDMD ***');
tic
A1 = X*Y'; A2 = Y*Y';
A_f = A1 * pinv(A2);
B1 = Y*X'; B2 = X*X';
A_b = B1 * pinv(B2);
A = (A_f/A_b)^0.5;
A = real(A);
toc

[Phi_sorted,D] = eig(A);
lambda = diag(D);
idx_exclude = (abs(angle(lambda)) < 2*pi*1.5*0.01) | (abs(lambda)>1);
lambda(idx_exclude) = [];
Phi_sorted(:,idx_exclude) = [];
[lambda,idx_sort] = sort(lambda,'descend');
Phi_sorted = Phi_sorted(:,idx_sort);
    

%%
save DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm Phi_sorted lambda A roi_exclude

%%
cd(current_path);
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');

atlasimage=niftiread('atlas/MMP_in_MNI_symmetrical_LR_diff.nii');
atlasinfo=niftiinfo('atlas/MMP_in_MNI_symmetrical_LR_diff.nii');

try
    label_idx_list;
catch
    label_idx_list = 1:N;
end

mkdir('DM_video_HCP_REST_fbDMD');
save_dir = [pwd filesep 'DM_video_HCP_REST_fbDMD'];

frame_dt = 0.5;

for pair_num = 1:9
    
    cd(save_dir);
    mkdir(['DM_pair', num2str(pair_num) '_4view']);
    cd(['DM_pair', num2str(pair_num) '_4view']);
    
    N=size(Phi_sorted,1);
    
    DM_conjugate1_num = 2*(pair_num-1)+1;
    DM_conjugate2_num = 2*pair_num;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    for frame = 0:1:100
        eigenstate = zeros(size(atlasimage));
        for roi=1:N 
            eigenstate(atlasimage==label_idx_list(roi))= real((lambda_conjugate1^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^(frame*frame_dt/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
        end
        
        atlasinfo.Datatype='double';
        niftiwrite(eigenstate,['DM_pair', num2str(pair_num), '_', num2str(frame,'%04d')],atlasinfo);
        
        fh = conn_mesh_display(['DM_pair', num2str(pair_num), '_', num2str(frame,'%04d'), '.nii']);
        fh('colormap','bluewhitered');
        fh('colorbar','on');
%         fh('colorbar','rescale',[-max(abs(Phi_sorted(:,DM_conjugate1_num)))/3,max(abs(Phi_sorted(:,DM_conjugate1_num)))/3]);
        fh('colorbar','rescale',[-max(abs(Phi_sorted(:,DM_conjugate1_num))),max(abs(Phi_sorted(:,DM_conjugate1_num)))]);
        fh('print',4,['DM_pair', num2str(pair_num), '_', num2str(frame,'%04d'),'.jpg'],'-r150','-nogui') 
        
        close;
    end
    
       
    % Define the folder where your images are stored
    imageFolder = pwd;
    imageFiles = dir(fullfile(imageFolder, '*.jpg')); % Adjust the pattern if needed
    
    niiFiles = dir(fullfile(imageFolder, '*.nii'));
    for ii = 1:length(niiFiles)
        delete([niiFiles(ii).folder,filesep,niiFiles(ii).name]);
    end

    % Create a VideoWriter object
    outputVideo = VideoWriter(fullfile(imageFolder, 'outputVideo.avi'));
    outputVideo.FrameRate = 10; % Set frame rate

    % Open the video writer
    open(outputVideo);

    % Loop through each image, read it, and write it to the video
    for i = 1:length(imageFiles)
        img = imread(fullfile(imageFolder, imageFiles(i).name));
        if i == 1
            desiredSize = [size(img,1), size(img,2)];
            resizedImg = img;
        else
            resizedImg = imresize(img, desiredSize); % Resize the image
        end
        
        % Calculate current time
        currentTime = 0 + (i - 1) * frame_dt;
        timeText = sprintf('t=%.1f', currentTime); % Format text

        % Add text to the image
        position = [20, 20]; % Position of the text (x, y)
        fontSize = 50; % Font size
        annotatedImg = insertText(double(resizedImg)/255, position, timeText, 'FontSize', fontSize, 'TextColor', 'black', 'BoxOpacity', 0);

        % Write the annotated image to the video
        writeVideo(outputVideo, annotatedImg);
    end

    % Close the video writer
    close(outputVideo);
end

%% individual fitting (based on optDMD)
cd(current_path);

num_DMs = 18;

U=Phi_sorted(:,1:num_DMs);

V=pinv(U);
D=zeros(num_DMs+1,length(tau));
for n=1:length(tau)
    if tau(n) == 0
        continue;
    end
    
    disp(['start: sub#' num2str(n)]);
    C_temp = zeros(num_DMs+1,num_DMs+1);
    B_temp = zeros(num_DMs+1,1);
    
    tic
    tau_n = sum(tau(1:n-1));
    X_temp = X(:,tau_n+1:tau_n+tau(n));
    Y_temp = Y(:,tau_n+1:tau_n+tau(n));
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
    toc

    D(:,n) = C_temp\B_temp;
    disp(['end: sub#' num2str(n)]);
end

save DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_indiv_18 tau Phi_sorted lambda D sub_ids num_DMs
