clear; clc;
current_path = pwd;
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered_meta.mat');
load('results/HCP_timeseries_cortical_subcortical_extracted_filtered.mat');

%% dFC - std
% MATLAB Script to Calculate Dynamic Functional Connectivity (dFC) Statistics: Std
% Author: Youngjo Song
% Date: 2024-10-11

%%% Define Parameters
tr = 0.72;
w = 45;
WL = round(w/tr);                  % Window length in TRs (adjust as needed)
num_subjects = size(time_series_denoised_filtered, 1); % 1096
num_sessions = size(time_series_denoised_filtered, 2); % 4
num_rois = 360;           % Number of ROIs
time_points = 1200;       % Number of time points per session

%%% Precompute all unique ROI pairs to avoid redundant computations
roi_pairs = nchoosek(1:num_rois, 2);
num_pairs = size(roi_pairs, 1);

% Initialize a matrix to store dstd for each subject, session, and ROI pair
% Dimensions: [Subjects x Sessions x ROI_Pairs]
dstd_all = zeros(num_subjects, num_sessions, num_pairs);

%%% Progress Tracking
fprintf('Calculating dFC Std for all subjects, sessions, and ROI pairs...\n');
total_iterations = num_subjects * num_sessions;
current_iteration = 0;

tic; % Start timer

%%% Loop Through Subjects and Sessions
for subj = 1:num_subjects
    for sess = 1:num_sessions
        current_iteration = current_iteration + 1;
        
        % Progress display every 100 iterations
        if mod(current_iteration, 100) == 0 || current_iteration == total_iterations
            elapsed_time = toc;
            est_total_time = (elapsed_time / current_iteration) * total_iterations;
            est_remaining = est_total_time - elapsed_time;
            fprintf('Processed %d/%d subjects (Elapsed Time: %.2f sec | Estimated Remaining: %.2f sec)\n', ...
                subj, num_subjects, elapsed_time, est_remaining);
        end
        
        % Check if the cell is not empty
        if ~isempty(time_series_denoised_filtered{subj, sess})
            % Extract fMRI data for the current subject and session
            % Data size: [ROIs x Time Points] = [360 x 1200]
            fmri_data = time_series_denoised_filtered{subj, sess}(1:num_rois, 1:time_points);
            
            % Transpose to [Time Points x ROIs] for consistency
            fmri_data = fmri_data';
            
            % high-pass filter
            fmri_data = highpass(fmri_data,1/(tr*WL),1/tr);
            
            % Handle Global Signal Regression if not already done
            % Uncomment the following lines if GSR needs to be performed here
            % group_data_gsr = regress_out_global_signal(fmri_data);
            
            % Loop through all unique ROI pairs
            parfor pair = 1:num_pairs
                roi1 = roi_pairs(pair, 1);
                roi2 = roi_pairs(pair, 2);
                
                % Extract time series for the ROI pair
                x = fmri_data(:, roi1);
                y = fmri_data(:, roi2);
                
                % Compute dstd using the defined function
                dstd = compute_dstd(x, y, WL);
                
                % Store the result
                dstd_all(subj, sess, pair) = dstd;
            end
        else
            % If the cell is empty, assign NaN to all ROI pairs
            dstd_all(subj, sess, :) = NaN;
        end
        disp(['sub #',num2str(subj,'%04d'),' finished']);
    end
end

dstd_group_vec = squeeze(nanmean(dstd_all,[1,2]));
dstd_group = zeros(num_rois);
for pair = 1:num_pairs
    roi1 = roi_pairs(pair, 1);
    roi2 = roi_pairs(pair, 2);
    
    dstd_group(roi1,roi2) = dstd_group_vec(pair);
    dstd_group(roi2,roi1) = dstd_group_vec(pair);
end

eval(['dsta_all_',num2str(w),' = dstd_all;']);
eval(['dstd_group_',num2str(w),' = dstd_group;']);

%%% Post-processing and Saving Results
save('results/dstd_all.mat', ['dstd_group_',num2str(w)],'-append');
disp('Dynamic Functional Connectivity (dFC) Std calculation completed successfully.');

%% dFC - Excursion
% MATLAB Script to Calculate Dynamic Functional Connectivity (dFC) Statistics: Excursion
% Author: Youngjo Song
% Date: 2024-10-12

%%% Define Parameters
WL = 30;                  % Window length in TRs (adjust as needed)
corr_threshold = 0.2;     % Correlation threshold for peak detection (optional)
num_subjects = size(time_series_denoised_filtered, 1); % 1096
num_sessions = size(time_series_denoised_filtered, 2); % 4
num_rois = 360;           % Number of ROIs
time_points = 1200;       % Number of time points per session

%%% Precompute all unique ROI pairs to avoid redundant computations
roi_pairs = nchoosek(1:num_rois, 2);
num_pairs = size(roi_pairs, 1);

% Initialize a matrix to store dstd for each subject, session, and ROI pair
% Dimensions: [Subjects x Sessions x ROI_Pairs]
execursion_all = zeros(num_subjects, num_sessions, num_pairs);

%%% Progress Tracking
fprintf('Calculating dFC Std for all subjects, sessions, and ROI pairs...\n');
total_iterations = num_subjects * num_sessions;
current_iteration = 0;

tic; % Start timer

%%% Loop Through Subjects and Sessions
for subj = 1:num_subjects
    for sess = 1:num_sessions
        current_iteration = current_iteration + 1;
        
        % Progress display every 100 iterations
        if mod(current_iteration, 100) == 0 || current_iteration == total_iterations
            elapsed_time = toc;
            est_total_time = (elapsed_time / current_iteration) * total_iterations;
            est_remaining = est_total_time - elapsed_time;
            fprintf('Processed %d/%d subjects (Elapsed Time: %.2f sec | Estimated Remaining: %.2f sec)\n', ...
                subj, num_subjects, elapsed_time, est_remaining);
        end
        
        % Check if the cell is not empty
        if ~isempty(time_series_denoised_filtered{subj, sess})
            % Extract fMRI data for the current subject and session
            % Data size: [ROIs x Time Points] = [360 x 1200]
            fmri_data = time_series_denoised_filtered{subj, sess}(1:num_rois, 1:time_points);
            
            % Transpose to [Time Points x ROIs] for consistency
            fmri_data = fmri_data';
            
            % Handle Global Signal Regression if not already done
            % Uncomment the following lines if GSR needs to be performed here
            % group_data_gsr = regress_out_global_signal(fmri_data);
            
            % Loop through all unique ROI pairs
            parfor pair = 1:num_pairs
                roi1 = roi_pairs(pair, 1);
                roi2 = roi_pairs(pair, 2);
                
                % Extract time series for the ROI pair
                x = fmri_data(:, roi1);
                y = fmri_data(:, roi2);
                
                % Compute Sliding Window Correlation
                [rho, ~] = compute_sliding_correlation(x, y, WL);

                % Compute Excursion
                excursion = compute_excursion(rho);
                
                % Store the result
                execursion_all(subj, sess, pair) = excursion;
            end
        else
            % If the cell is empty, assign NaN to all ROI pairs
            execursion_all(subj, sess, :) = NaN;
        end
        disp(['sub #',num2str(subj,'%04d'),' finished']);
    end
end

execursion_group_vec = squeeze(nanmean(execursion_all,[1,2]));
execursion_group = zeros(num_rois);
for pair = 1:num_pairs
    roi1 = roi_pairs(pair, 1);
    roi2 = roi_pairs(pair, 2);
    
    execursion_group(roi1,roi2) = execursion_group_vec(pair);
    execursion_group(roi2,roi1) = execursion_group_vec(pair);
end

%%% Post-processing and Saving Results
save('execursion_all.mat', 'execursion_group','execursion_all');
disp('Dynamic Functional Connectivity (dFC) Execursion calculation completed successfully.');

%% functions

function excursion = compute_excursion(rho)
% compute_excursion Calculates the Excursion of dynamic functional connectivity (dFC Excursion)
%
% Syntax:
%   excursion = compute_excursion(rho)
%
% Inputs:
%   rho - A vector of sliding window correlation coefficients
%
% Outputs:
%   excursion - The Excursion measure as defined by Zalesky et al. (2014)
%
% Example:
%   rho = randn(100,1); % Example correlation coefficients
%   excursion = compute_excursion(rho);

    % Compute the median of the correlation coefficients
    median_rho = median(rho);
    
    % Determine the sign of (rho - median)
    sign_rho = sign(rho - median_rho);
    
    % Find median crossing points where the sign changes
    sign_diff = diff(sign_rho);
    crossing_points = find(sign_diff ~= 0) + 1; % +1 to get the index after the sign change
    
    % Include the start and end points as crossing points
    crossing_points = [1; crossing_points; length(rho) + 1]; % length(rho)+1 to include the last point
    
    % Initialize Excursion
    excursion = 0;
    
    % Define alpha and beta
    alpha = 0.9;
    beta = 1;
    
    % Iterate through each excursion
    for k = 1:length(crossing_points)-1
        tk = crossing_points(k);
        tk1 = crossing_points(k+1);
        
        % Excursion length
        lk = tk1 - tk;
        
        % Excursion height (maximum rho in the excursion)
        if tk1 > length(rho)
            rho_segment = rho(tk:end);
        else
            rho_segment = rho(tk:tk1-1);
        end
        hk = max(rho_segment);
        
        % Accumulate the weighted excursion
        excursion = excursion + (lk^alpha) * (hk^beta);
    end
end

function [rho, rho_mean] = compute_sliding_correlation(x, y, WL)
% compute_sliding_correlation Computes sliding window Pearson correlation between two time series
%
% Syntax:
%   [rho, rho_mean] = compute_sliding_correlation(x, y, WL)
%
% Inputs:
%   x  - A vector of length N representing the time series of ROI 1
%   y  - A vector of length N representing the time series of ROI 2
%   WL - Window length (number of consecutive time points)
%
% Outputs:
%   rho      - Vector of Pearson correlation coefficients for each window
%   rho_mean - Mean of the correlation coefficients
%
% Example:
%   N = 1200;
%   x = randn(N, 1);
%   y = randn(N, 1);
%   WL = 30;
%   [rho, rho_mean] = compute_sliding_correlation(x, y, WL);

    % Ensure inputs are column vectors
    x = x(:);
    y = y(:);
    
    % Check that x and y are of the same length
    if length(x) ~= length(y)
        error('Time series x and y must be of the same length.');
    end
    
    N = length(x);
    
    % Number of sliding windows
    num_windows = N - WL + 1;
    
    % Initialize vector to store correlation coefficients
    rho = zeros(num_windows, 1);
    
    % Vectorized computation for efficiency
    % Create matrices where each row is a window for x and y
    X = zeros(num_windows, WL);
    Y = zeros(num_windows, WL);
    for w = 1:WL
        X(:, w) = x(w:(w + num_windows -1));
        Y(:, w) = y(w:(w + num_windows -1));
    end
    
    % Compute mean for each window
    mean_X = mean(X, 2);
    mean_Y = mean(Y, 2);
    
    % Demean the data
    X_demean = X - mean_X;
    Y_demean = Y - mean_Y;
    
    % Compute numerator and denominator for Pearson correlation
    numerator = sum(X_demean .* Y_demean, 2);
    denominator = sqrt(sum(X_demean.^2, 2) .* sum(Y_demean.^2, 2));
    
    % Avoid division by zero
    denominator(denominator == 0) = eps;
    
    % Compute Pearson correlation coefficients
    rho = numerator ./ denominator;
    
    % Compute mean correlation across all windows
    rho_mean = mean(rho);
end