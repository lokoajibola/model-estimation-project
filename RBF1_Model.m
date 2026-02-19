function rbf1 = RBF1_Model(Nx)
%     Nx2 = Nx;
    Nx = 2;
    %Mahalanobis distance equation
    rbf1.Mahal_dist = @(a1,m,c)sqrt((a1-m)*inv(c)*(a1-m)');
%     rbf1.Mahal_dist2 = @(a1,m,c)sqrt((a1-m)*inv(c)*(a1-m)');
    %initialize variables as vectors and matrix in symbol form
    a1 = sym('a1', [1 Nx]);
%     m1 = sym('m1', [1 Nx1]);
    m = sym('m', [1 Nx]);
    c = sym('c', [Nx Nx]);
    
    %compute the derivatives
    rbf1.Jac = jacobian(sqrt((a1-m)*inv(c)*(a1-m)'),m);
    rbf1.Hes = hessian(sqrt((a1-m)*inv(c)*(a1-m)'),m);
    
%     rbf1.Jac2a = jacobian(rbf1.Jac(1),m(1));
%     if Nx == 2
%         rbf1.Jac2b = jacobian(rbf1.Jac(2),m(2));
%     end
    
    %RBF functions for working with options prices
    rbf1.RBF = @(real_input,X_p,cov,wts)get_C_X(real_input, wts, X_p, cov);    
    rbf1.choice = 'md_RBF';
end