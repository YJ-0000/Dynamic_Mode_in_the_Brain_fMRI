function [F,p_val] = F_test_GLM(X,y,C)
    J = size(X,1);
    p = rank(X);
    
    C0 = eye(size(X,2)) - C*pinv(C);
    
    X0 = X * C0;
    
    p2 = rank(X0);
    
    p1 = p - p2;
    
    
    R0 = eye(J) - X0 * pinv(X0);
    R = eye(J) - X * pinv(X);
    
    M = R0 - R;
    
    F = ((J-p)/p1) * (y' * M * y) / (y' * R * y);
    p_val = 1 - fcdf(F, p1, J-p);
end