% Generate synthetic data to analyse the determinant of the hessian matrix

clear
close all


Data.Model_choice = 'BS';
Data.Filter_choice = 'SOEKF';

Data.Data_type = 'real';
% Data.Data_type = 'synthetic';

Data.data_yr = 2014;    % 2014 - 6800 or 6700; 2013 - 6100 or 6000; 2012 - 5500 or 5600
Data.X = 6700;

Data.option_type = 'call';
% Data.option_type = 'put';

% get model
Model = Black_Scholes_Model(2);
Data = get_data(Model, Data);

% generate S/X and T-t using mesh grid for the ratio values
t_m = (linspace(1,0,length(Data.real_in)))';
S_X  = (linspace(0.8,1.2,length(Data.real_in)))';
X = 6600;
S = S_X.*X;
[t_m_mesh, S_X_mesh] = meshgrid(t_m,S_X);

if strcmp(Data.Data_type,'real')
    X = Data.X;
    S = Data.S;
    t_m = Data.t_m;
    S_X = S./X;
    [t_m_mesh, S_X_mesh] = meshgrid(t_m,S_X);
end

r = Data.r;
sig = Data.sig;

for i = 1:200 %length(t_m_mesh)-1 %
    for j = 1:200 %length(S_X)-1 %
        Hes = Model.Hes1(S(j), X, r(j), sig(j), t_m(i));
        Jac = Model.Jac2(S(j), X, r(j), sig(j), t_m(i));
        det_Hes(j,i) = det((Hes));
        det_Jac(j,i) = det(Jac);
    end
end


det_Hes2 = det_Hes./(1+(det_Jac.^2)).^1.5;

figure
imagesc(det_Hes)
