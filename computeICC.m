function ICC = computeICC(y1, y2)
    % computeICC calculates the Intraclass Correlation Coefficient (ICC)
    % using a two-way mixed-effects model (ICC(2,1) Absolute Agreement).
    %
    % Inputs:
    %   y1 - Vector of measurements from session 1
    %   y2 - Vector of measurements from session 2
    %
    % Output:
    %   ICC - Intraclass Correlation Coefficient value

    if length(y1) ~= length(y2)
        error('Input vectors must be of the same length.');
    end

    % Combine data into a matrix
    data = [y1(:), y2(:)]; % Ensure column vectors

    nTargets = size(data, 1); % Number of targets (subjects)
    nRatters = size(data, 2); % Number of raters (sessions)

    % Mean per target
    targetMeans = mean(data, 2);
    % Mean per rater
    raterMeans = mean(data, 1);
    % Grand mean
    grandMean = mean(data(:));

    % Sum of Squares
    SStotal = sum((data(:) - grandMean).^2);
    SSrater = nTargets * sum((raterMeans - grandMean).^2);
    SStarget = nRatters * sum((targetMeans - grandMean).^2);
    SSerror = SStotal - SSrater - SStarget;

    % Degrees of freedom
    dfTotal = nTargets * nRatters - 1;
    dfRater = nRatters - 1;
    dfTarget = nTargets - 1;
    dfError = dfTarget * dfRater;

    % Mean Squares
    MSrater = SSrater / dfRater;
    MStarget = SStarget / dfTarget;
    MSerror = SSerror / dfError;

    % ICC(2,1) Absolute Agreement
    ICC = (MStarget - MSerror) / (MStarget + (nRatters - 1) * MSerror + (nRatters * (MSrater - MSerror)) / nTargets);
end