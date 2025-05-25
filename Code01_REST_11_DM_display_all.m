clear;clc;
current_path = pwd;
conn;
close all;
load DMs/DM_cortical_subcortical_ext_fbDMD_noROInorm_subExclude.mat
%%
cd(current_path);
labels = cifti_read('atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
labels_subcortical = cifti_read('atlas/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii');
labels_subcortical_data = labels_subcortical.cdata;

atlasimage=niftiread('atlas/MMP_in_MNI_symmetrical_LR_diff.nii');
atlasinfo=niftiinfo('atlas/MMP_in_MNI_symmetrical_LR_diff.nii');

N_cortex = 360;
atlasimage_temp = atlasimage;
for n_roi = 1:(N_cortex/2)
    atlasimage(atlasimage_temp==n_roi) = n_roi + 180;
end
for n_roi = ((N_cortex/2)+1):N_cortex
    atlasimage(atlasimage_temp==n_roi) = n_roi - 180;
end

label_idx_list_cortex = 1:N_cortex;

label_idx_list_subcortical = unique(labels_subcortical_data(:));
label_idx_list_subcortical(label_idx_list_subcortical==0) = [];

label_idx_list_subcortical(roi_exclude) = [];

mkdir('DM_video_HCP_REST_fbDMD');
save_dir = [pwd filesep 'DM_video_HCP_REST_fbDMD'];

frame_dt = 0.5;
TRtarget = 1.5;

roi_list = {'cortex','hippocampus','amygdala','thalamus','striatum','brainstem','cerebellum'};
roi_names_display = {'Cortex','Hippocampus','Amygdala','Thalamus','Striatum','Brainstem','Cerebellum'};


%% colormap for subcortex and cerebellum
% Define the number of points for each segment
n = 128; % Half of 256

% -------------------------------
% Define Red to Grey Transition
% -------------------------------

% Red channel: 1 to 0.5
R_red_to_grey = linspace(1, 0.5, n);

% Green channel: 0 to 0.5
G_red_to_grey = linspace(0, 0.5, n);

% Blue channel: 0 to 0.5
B_red_to_grey = linspace(0, 0.5, n);

% -------------------------------
% Define Grey to Blue Transition
% -------------------------------

% Red channel: 0.5 to 0
R_grey_to_blue = linspace(0.5, 0, n);

% Green channel: 0.5 to 0
G_grey_to_blue = linspace(0.5, 0, n);

% Blue channel: 0.5 to 1
B_grey_to_blue = linspace(0.5, 1, n);

% -------------------------------
% Combine the Segments
% -------------------------------

% Concatenate the segments for each RGB channel
R = [R_red_to_grey, R_grey_to_blue];
G = [G_red_to_grey, G_grey_to_blue];
B = [B_red_to_grey, B_grey_to_blue];

% Combine into a single colormap matrix
cmap = [R', G', B'];
cmap = flipud(cmap);

% -------------------------------
% Apply Scaling Transformation (Optional)
% -------------------------------

% The original code applied a scaling transformation to enhance the gradient.
% We'll apply a similar transformation adjusted for 256 points.

% % Define the scaling factor across the range [-1, 1]
% scale_factor = abs(linspace(-1, 1, 256))';
% 
% % Apply the transformation
% cmap = scale_factor .* cmap + (1 - scale_factor) * 1;

% cmap = slanCM(1, 256);
% Scales
%%
for pair_num = 1:12
    
    if pair_num == 1
        ref_t = 27;
    elseif pair_num == 3
        ref_t = 31+30;
    elseif pair_num == 5
        ref_t = 43;
    elseif pair_num == 2
        ref_t = 56.5 + 19;
    elseif pair_num == 4
        ref_t = 14 + 37.5;
    else
        ref_t = 0;
    end
    
    cd(save_dir);
    mkdir(['DM_pair', num2str(pair_num), '_all']);
    cd(['DM_pair', num2str(pair_num), '_all']);
    
    target_dir = pwd;
    
    N=size(Phi_sorted,1);
    
    DM_conjugate1_num = 2*(pair_num-1)+1;
    DM_conjugate2_num = 2*pair_num;
    
    lambda_conjugate1 = lambda(DM_conjugate1_num);
    lambda_conjugate2 = lambda(DM_conjugate2_num);
    
    if pair_num == 1
        min_scale = -max(abs(Phi_sorted(:,DM_conjugate1_num)))/3;
        max_scale = max(abs(Phi_sorted(:,DM_conjugate1_num)))/3;
    else
        min_scale = -max(abs(Phi_sorted(:,DM_conjugate1_num)));
        max_scale = max(abs(Phi_sorted(:,DM_conjugate1_num)));
    end
    
    frame_num = ceil(TRtarget * 2*pi / abs(angle(lambda_conjugate1))  / frame_dt);
    
    for frame = 0:1:frame_num
        for n_roi = 1:length(roi_list)
            target_roi = roi_list{n_roi};
            
            if strcmp(target_roi,'cortex')
                eigenstate = zeros(size(atlasimage));
                for roi=1:length(label_idx_list_cortex)
                    eigenstate(atlasimage==label_idx_list_cortex(roi)) ...
                        = real(((lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget))/abs(lambda_conjugate1^(ref_t/TRtarget)))*Phi_sorted(roi,DM_conjugate1_num) ...
                            + ((lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget))/abs(lambda_conjugate2^(ref_t/TRtarget)))*Phi_sorted(roi,DM_conjugate2_num));
                end
                atlasinfo.Datatype='double';
                niftiwrite(eigenstate,'temp',atlasinfo);

                fh = conn_mesh_display('temp.nii');
                fh('colormap','bluewhitered');
                fh('background',[1,1,1]);
                fh('colorbar','on');
                fh('colorbar','rescale',[min_scale,max_scale]);
                fh('print',4,['DM_',target_roi,'_pair', num2str(pair_num), '_', num2str(frame,'%04d'),'.jpg'],'-r150','-nogui') 

                close;
            else
                eigenstate_cifti = zeros(size(labels_subcortical.cdata));
                for roi=1:length(label_idx_list_subcortical)
                    eigenstate_cifti(labels_subcortical_data==label_idx_list_subcortical(roi)) = ...
                        real((lambda_conjugate1^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate1_num) + (lambda_conjugate2^((frame*frame_dt+ref_t)/TRtarget))*Phi_sorted(roi,DM_conjugate2_num));
                end
                cd(current_path);
                plot_subcortex(eigenstate_cifti,[target_dir,filesep,'DM_',target_roi,'_pair', num2str(pair_num), '_', num2str(frame,'%04d'),'.jpg'], ...
                            cmap,  min_scale, max_scale, target_roi);
                close;
            end
            
        end
        
        %%
        positions = [
            0.00, 0.3, 0.7, 0.7;  % Position for image 1
            0.00, 0.00, 0.25, 0.3;  % Position for image 2
            0.22, 0.00, 0.25, 0.3;  % Position for image 3
            0.45, 0.00, 0.25, 0.3;  % Position for image 4
            0.7, 0.7, 0.25, 0.3;  % Position for image 5
            0.7, 0.40, 0.25, 0.3;  % Position for image 6
            0.6, 0, 0.45, 0.45;  % Position for image 7
        ];
        title_positions = [
            0.43, 0.9;
            0.5, 0.9;
            0.5, 0.9;
            0.5, 0.9;
            0.5, 0.9;
            0.5, 0.9;
            0.5, 0.81;
        ];
        
        cd(target_dir);
        figure('Position', [100, 100, 900, 700]);
        % Loop through the seven images
        for n_roi = [7,1:6]
            target_roi = roi_list{n_roi};
            % Construct the filename for each image (assuming names as 'image1.jpg' to 'image7.jpg')
            filename = ['DM_',target_roi,'_pair', num2str(pair_num), '_', num2str(frame,'%04d'),'.jpg'];

            % Load the image
            img = imread(filename);

            % Create axes with specified position
            ax = axes('InnerPosition', positions(n_roi, :));

            % Display the image within the axes
            imshow(img, 'Parent', ax);

            % Optionally, turn off the axis for the image
            axis(ax, 'off');
            
            text(title_positions(n_roi,1), title_positions(n_roi,2), target_roi, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'Parent', ax);
        end
        saveas(gcf, ['DM_all_pair', num2str(pair_num), '_', num2str(frame,'%04d'),'.jpg']);
%%
        close;
        
    end
    
       
    % Define the folder where your images are stored
    imageFolder = target_dir;
    imageFiles = dir(fullfile(imageFolder, 'DM_all*.jpg')); % Adjust the pattern if needed
    
    niiFiles = dir(fullfile(imageFolder, '*.nii'));
    for ii = 1:length(niiFiles)
        delete([niiFiles(ii).folder,filesep,niiFiles(ii).name]);
    end
    for n_roi = 1:7
        redunFiles = dir(fullfile(imageFolder, ['DM_', roi_list{n_roi},'*.jpg']));
        for ii = 1:length(redunFiles)
            delete([redunFiles(ii).folder,filesep,redunFiles(ii).name]);
        end
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
        timeText = sprintf('t=%.1fs', currentTime); % Format text

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