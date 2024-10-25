function [t,p_val] = t_test_GLM(X,y,c,d)
    if nargin < 4
        d = 0;
    end
    
    J = size(X,1);
    p = rank(X);
    
    beta = pinv(X) * y;
    
    resid_e = y - X * beta;
    
    sigma_squared = (resid_e' * resid_e)/(J-p);
    
    t = (c'*beta - d)/sqrt(sigma_squared * c' * pinv(X'*X) * c);
    p_val = 2 * (1 - tcdf(abs(t), J-p));
end