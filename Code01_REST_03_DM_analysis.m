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
    disp(ii);
    nsub = ceil(ii/4); nses = rem(ii,4); if nses==0; nses=4;end
    if isempty(time_series_denoised_filtered{nsub,nses})
        continue
    end
    if isnan(sum(time_series_denoised_filtered{nsub,nses},'all'))
        continue
    end
    i_num = i_num + 1;
end

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
        if ~isempty(time_series_denoised_filtered{nsub,nses}) && ~isnan(sum(time_series_denoised_filtered{nsub,nses},'all'))
            y = time_series_denoised_filtered{nsub,nses};
            if t_sample ~= TRtarget
                pp = spline(t, y);          % Compute piecewise polynomial (B-spline) representation
                y_fine = ppval(pp, t_fine); % Evaluate the piecewise polynomial at the finer time points
            else
                y_fine = y;
            end
            
            % ROI normalization
%             for nroi = 1:size(y_fine,1)
%                 y_fine(nroi,:) = normalize(y_fine(nroi,:));
%             end
            
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
% X(:,i_num*(length(t_fine)-1)+1:end) = [];
% Y(:,i_num*(length(t_fine)-1)+1:end) = [];


clear time_series_denoised_filtered

%%
A=X/Y;
B=Y/X;
A = (A/B)^0.5;
A = real(A);

[Phi,Diag]=eig(A);
dd = abs(diag(Diag));
[dd_sorted,idx_sorted] = sort(dd,'descend');
max_number_eigenstates=knee_pt(dd_sorted);
max_number_eigenstates=2*floor(max_number_eigenstates/2);
Phi_sorted=Phi(:,idx_sorted);

lambda = diag(Diag);
lambda = lambda(idx_sorted);
D_ang = angle(lambda);

aaa = log(0.5) ./ log(abs(lambda));
bbb = 2*pi./D_ang;
ccc = (aaa./abs(bbb));

%%
save DMs/DM_cortical_subcortical Phi_sorted lambda A max_number_eigenstates

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

mkdir('DM_video_HCP_REST');
save_dir = [pwd filesep 'DM_video_HCP_REST'];

frame_dt = 0.5;

for pair_num = 1:1%round(max_number_eigenstates/2)
    
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
%         fh('colorbar','rescale',[-max(abs(Phi_sorted(:,DM_conjugate1_num)))/2,max(abs(Phi_sorted(:,DM_conjugate1_num)))/2]);
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

%% individual fitting

U=Phi_sorted(:,(1:max_number_eigenstates));
V=pinv(U);

% C=zeros(max_number_eigenstates,max_number_eigenstates,length(tau));
% B=zeros(max_number_eigenstates,length(tau));
D=zeros(max_number_eigenstates,length(tau));
for n=1:length(tau)
    if tau(n) == 0
        continue;
    end
    
    disp(['start: sub#' num2str(n)]);
    C_temp = zeros(max_number_eigenstates,max_number_eigenstates);
    B_temp = zeros(max_number_eigenstates,1);
    
    tic
    tau_n = sum(tau(1:n-1));
    X_temp = X(:,tau_n+1:tau_n+tau(n));
    Y_temp = Y(:,tau_n+1:tau_n+tau(n));
    VY = V * Y_temp;
    for i=1:max_number_eigenstates
        temp1 = U(:,i) * VY(i,:);
        for j=1:max_number_eigenstates  
            temp2 = U(:,j) * VY(j,:);
            C_temp(i,j) = 2 * sum(dot(temp1, temp2, 1));
        end
    end
    for k=1:max_number_eigenstates
        B_temp(k) = sum(diag(2*X_temp'*U(:,k)*VY(k,:)));
    end
    toc

%     C(:,:,n) = C_temp;
%     B(:,n) = B_temp;
    D(:,n) = C_temp\B_temp;
    disp(['end: sub#' num2str(n)]);
end

save results/DM_cortical_subcortical_indiv tau Phi_sorted lambda A D sub_ids max_number_eigenstates

