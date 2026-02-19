function [Zh, Jac, Hes, X_p] = eval_BS(input, Model)
% This function uses the input parameters to evaluate:
% the Black-Scholes Call option price estimate [Zh]
% the Jacobian of the Black-Scholes [Jac]
% the Hessian of the Black-Scholes [Hes]

k = input.k - 1;
S = input.S(k);
X = input.X;    % Strike price 
t_m = input.real_in(k,2);
X_p = input.X_a;    % State vector
[Nx1, Nx2] = size(X_p);

if Nx2 == 1
    sig = X;
    r = input.r(k);
    Jac = Model.dCs(S,X,r,sig,t_m); % evaluate the Jacobian
    Hes = Model.Hes1(S,X,r,sig,t_m); % evaluate the Hessian
    
elseif Nx2 == 2
    sig = X_p(1); r = X_p(2);
    
    syms sig r;
    Jac_eq = jacobian(Model.C, [sig r]);
    Hes_eq = hessian(Model.C, [sig r]);
    sig1 = X_p(1); r1 = X_p(2);
    Jac = eval(subs(Jac_eq, [S X t_m sig r],[S X t_m sig1 r1])); % evaluate the Jacobian
    Hes = eval(subs(Hes_eq, [S X t_m sig r],[S X t_m sig1 r1])); % evaluate the Hessian
    
    
%     Jac = Model.Jac1(S,X,r,sig,t_m); 
%     Hes = Model.Hes1(S,X,r,sig,t_m); 
%     %Hes2 = Model.Hes2(S,X,r,sig,t_m);
end

% measurement estimate
Zh = eval(subs(Model.C, [S X t_m sig r],[S X t_m sig1 r1]));

end