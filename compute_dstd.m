function dstd = compute_dstd(x, y, WL)
% compute_dstd Calculates the standard deviation of dynamic functional connectivity (dFC Std)
%
% Syntax:
%   dstd = compute_dstd(x, y, WL)
%
% Inputs:
%   x  - A vector of length N representing the time series of ROI 1
%   y  - A vector of length N representing the time series of ROI 2
%   WL - Window length (number of consecutive time points)
%
% Outputs:
%   dstd - The standard deviation of the sliding window correlation coefficients
%
% Example:
%   N = 1200;
%   x = randn(N, 1);
%   y = randn(N, 1);
%   WL = 30;
%   dstd = compute_dstd(x, y, WL);

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
    
    % Compute mean and standard deviation for x and y within each window
    % Vectorized computation for efficiency
    % Create a matrix where each row is a window for x and y
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
    
    % Compute standard deviation as per the formula
    dstd = sqrt(mean((rho - rho_mean).^2));
    
end