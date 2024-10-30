function R2_accum = linear_LASSO_yj_training(Y_fitting,X_fitting,Y_validation,X_validation,lambda_list)
    c0 = sum(dot(Y_fitting, Y_fitting, 1));
    b0 = sum(dot(Y_fitting, X_fitting, 1));
    a_null = real(b0/c0);
    
    s_m = X_fitting - a_null*Y_fitting;
    
    Theta = zeros(size(Y_fitting, 1), size(Y_fitting, 1),length(lambda_list));
    parfor i = 1:size(Y_fitting, 1)
        % Transpose data to match lasso's expected input
        [B, ~] = lasso(Y_fitting', s_m(i,:)', 'lambda', lambda_list, 'Standardize', true);
        % If multiple coefficients are returned, select the one corresponding to lambda
        Theta(i, :, :) = B;  % Assuming lambda is a scalar
    end
    
    R2_accum = zeros(length(lambda_list),1);
    for n_inner = 1:length(lambda_list)
        W = (a_null * eye(size(Theta,1)) + squeeze(Theta(:,:,n_inner)));
        Y_test_plus_hat = W  * Y_validation;
        E_test = X_validation - Y_test_plus_hat;                                     % Prediction error
        R2 = mean(1 - sum(E_test.^2, 2) ./ sum((X_validation - mean(X_validation, 2)).^2, 2));
        R2_accum(n_inner) = R2;
    end
end

