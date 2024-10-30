function R2 = linear_LASSO_yj(Y_fitting,X_fitting,Y_validation,X_validation,Lambda)
    c0 = sum(dot(Y_fitting, Y_fitting, 1));
    b0 = sum(dot(Y_fitting, X_fitting, 1));
    a_null = real(b0/c0);
    
    s_m = X_fitting - a_null*Y_fitting;
    
    Theta = zeros(size(Y_fitting, 1));
    parfor i = 1:size(Y_fitting, 1)
        % Transpose data to match lasso's expected input
        [B, ~] = lasso(Y_fitting', s_m(i,:)', 'lambda', Lambda, 'Standardize', true);
        % If multiple coefficients are returned, select the one corresponding to lambda
        Theta(i, :) = B;  % Assuming lambda is a scalar
    end
    W = (a_null*eye(size(Theta)) + Theta) ;
    Y_test_plus_hat = W  * Y_validation;
    E_test = X_validation - Y_test_plus_hat;                                     % Prediction error
    R2 = mean(1 - sum(E_test.^2, 2) ./ sum((X_validation - mean(X_validation, 2)).^2, 2));
end

