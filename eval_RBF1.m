function [Zh, Jac, Hes] = eval_RBF1(input, Model)
% This function uses the input parametrs in a RBF Model to evaluate:
% the  Call option price estimate [Zh]
% the Jacobian [Jac]
% the Hessian  [Hes]

k = input.k-1;
c1 = input.cov;
if k == 2
    X_p = input.X_a;
else
    X_p = (input.X_a);
end
[Nx1, Nx2] = size(X_p);
real_input = input.model_in(k,:);

% Evaluate the Jacobian
Jac = zeros(Nx1,Nx2);
for count = 1:Nx1
    c2 = input.cov;
    Jac(count,:) = eval_jac(real_input,X_p(count,:),c2,Model.Jac);
    % Evaluate the Hessian
    Hes{count,:} = eval_jac(real_input,X_p(count,:),c2,Model.Hes);
end


% Measurement estimation
if strcmp(input.Data_type, 'real')
    Zh = Model.RBF(real_input, X_p, c1, input.RBFwts);
elseif strcmp(input.Data_type, 'synthetic')
    seq = 1:size(X_p,1);
    lamda = seq./(sum(seq));
    for i = 1:size(X_p,1)
        Zh1(i,:) = lamda(i) * Model.Mahal_dist(real_input, X_p(i,:), c1);
    end
    Zh = sum(Zh1);
end
Jac = Jac';
end