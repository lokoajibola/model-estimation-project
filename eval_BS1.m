function [Zh, Jac, Hes, X_p] = eval_BS1(input, Model)
% This function uses the input parameters to evaluate:
% the Black-Scholes Call option price estimate [Zh]
% the Jacobian of the Black-Scholes [Jac]
% the Hessian of the Black-Scholes [Hes]

k = input.k-1;
S = input.S(k);
X = input.X;    % Strike price
t_m = input.real_in(k,2);
X_p = input.X_a;    % State vector
[~ , Nx2] = size(X_p);

if Nx2 == 1
    sig = X_p;
    r = input.r(k);
    Jac = Model.dCs(S,X,r,sig,t_m);             % evaluate the Jacobian
    Hes = Model.Hes1(S,X,r,sig,t_m);            % evaluate the Hessian
    if strcmp(input.option_type,'put')
        Jac = Model.dPs(S,X,r,sig,t_m);         % evaluate the Jacobian
        Hes = Model.Hes2(S,X,r,sig,t_m);        % evaluate the Hessian
    end
    
elseif Nx2 == 2
    sig = (X_p(1)); r = (X_p(2));
    if strcmp(input.Data_type,'real')
        sig = exp(X_p(1)); r = exp(X_p(2));
    end
    Jac = Model.Jac1(S,X,r,sig,t_m);            % evaluate the Jacobian
    if strcmp(input.option_type,'put')
        Jac = Model.Jac3(S,X,r,sig,t_m);        % evaluate the Jacobian
    end
    Hes{1} = Model.Hes1(S,X,r,sig,t_m);         % evaluate the Hessian
    Hes{2} = Model.Hes2(S,X,r,sig,t_m);
    
end

% measurement estimate
Zh = Model.C(S,X,r,sig,t_m);
if strcmp(input.option_type,'put')
    Zh = Model.P(S,X,r,sig,t_m);
end

end