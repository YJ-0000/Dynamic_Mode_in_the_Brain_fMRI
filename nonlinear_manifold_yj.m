function R2 = nonlinear_manifold_yj(Y_fitting,X_fitting,Y_validation,X_validation,h)
    K_h = @(d, h) exp(-d.^2/(2*(h^2)));  
    
    h_orig = h;
    base_med_dist = 50;
    
    c0 = sum(dot(Y_fitting, Y_fitting, 1));
    b0 = sum(dot(Y_fitting, X_fitting, 1));
    a_null = real(b0/c0);
    
    s_m = X_fitting - a_null*Y_fitting;
    

    Phi = [ones(1, size(Y_fitting,2), size(Y_validation,2)); Y_fitting - permute(Y_validation, [1 3 2])]; % Matrix of regressors
    dist = permute(sqrt(sum(Phi(2:end, :, :).^2, 1)), [3 2 1]); 
    
    h = h_orig * median(dist(:)) / base_med_dist;
    
    K = K_h(dist, h);
    
    Theta0 = nan(size(Y_validation));                                                    % The matrix of parameters. Each slice is for one test point. Within each slice, the first column is the zero'th order Taylor term (at the corresponding test point) and the remaining n x n matrix is the coefficient of the linear Taylor term (which is discarded).
    parfor i_test = 1:size(Y_validation,2)
        G = Phi(:, :, i_test) * diag(K(i_test, :)) * Phi(:, :, i_test)';
        g = s_m * diag(K(i_test, :)) * Phi(:, :, i_test)';
        theta = g * pinv(G);
        Theta0(:, i_test) = theta(:, 1);
    end
    Y_test_plus_hat = a_null*Y_validation + Theta0;    
    E_test = X_validation - Y_test_plus_hat;                                     % Prediction error
    R2 = mean(1 - sum(E_test.^2, 2) ./ sum((X_validation - mean(X_validation, 2)).^2, 2));
end

