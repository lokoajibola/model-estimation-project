function [input] = KF(Model,input, k)

input.k = k;
Q = input.Q;
R = input.R;
X_a = input.X_a;
P_a = input.P_a;
Z = input.model_out(k-1);

[Nx1, Nx2] = size(X_a);

%% PREDICTION
X_p = X_a;      % state prediction
P_p = P_a + Q;  % prior covariance matrix

switch Model.choice
    case 'BS'
        [Zh, Jac, Hes, X_p] = eval_BS(input,Model);  % Evaluation of the derivatives
        S1 = Jac*P_p*Jac' + R;
        
        if strcmp (input.Filter_choice,'SOEKF')
            Hes1 = Hes{1};
            Hes2 = Hes{2};
            S1 = Jac*P_p*Jac' + 0.5*trace(Hes1*P_p) + R;
            input.det_Hes(k,:) = det(Hes1);
            input.det_Hes2(k,:) = max(eig(Hes1));
            
            if strcmp(input.option_type,'put')
                S1 = Jac*P_p*Jac' + 0.5*trace(Hes2*P_p) + R;
                input.det_Hes(k,:) = det(Hes2);
                input.det_Hes2(k,:) = max(eig(Hes2));
            end
            input.det_S1_ekf(k,:) = det(Jac*P_p*Jac');
            
        end
        
    case 'md_RBF'
        [Zh, Jac, Hes] = eval_RBF1(input, Model);  %% Evaluation of the derivatives
        S1 = Jac*P_p*Jac' + R;
        
        if strcmp (input.Filter_choice,'SOEKF')
            for i = 1:size(X_a,1)
                Hes1(i,i) = Hes{i,:}(1,1);
                Hes2(i,i) = Hes{i,:}(2,2);
            end
            new_Hes1 = [trace(Hes1*P_p*Hes1*P_p),trace(Hes1*P_p*Hes2*P_p);trace(Hes2*P_p*Hes1*P_p),trace(Hes2*P_p*Hes2*P_p)] ;
            S1 = Jac*P_p*Jac' + 0.5*(new_Hes1) + R;
            input.det_Hes(k,:) = det(new_Hes1);
            input.det_Hes2(k,:) = eig(new_Hes1)';
            input.det_S1_ekf(k,:) = det(Jac*P_p*Jac');
        end
end


K = P_p*Jac'*inv(S1);   % kalman gain
v = (Z - Zh);           % measurement residual

%% CORRECTION
switch Model.choice
    case 'BS'
        % state update
        input.X_a = X_p + K'*v;
        
        % covariance update
        input.P_a = (eye(Nx2) - K*Jac)*P_p';
        if strcmp (input.Filter_choice,'SOEKF')
            input.P_a = (eye(Nx2) - K*Jac)*P_p*(eye(Nx2) - K*Jac)' + (K*R*K');
        end
        
    case 'md_RBF'
        % state update
        input.X_a = X_p + K*v;
        
        %covariance update
        input.P_a = (eye(Nx1) - K*Jac)*P_p;
        if strcmp (input.Filter_choice,'SOEKF')
            input.P_a = (eye(Nx1) - K*Jac)*P_p*(eye(Nx1) - K*Jac)' + (K*R*K');
        end
        
        input.X_p = X_p;
        input.Jac{k,:} = Jac;
end

end