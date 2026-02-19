clear, close all

% This is the main program. This code illustrates the use of Extended Kalman
% filter (Both first and second order) on the estimation of call and Put option
% prices. This is done by
% [1] Estimating the volatility and the interest rates as state parameters
% for the Black-Scholes model
% [2] Estimating the mean as state parameters for the RBF model

% The real data uses FTSE100 prices, FTSE100 call prices and UK
% gilts as risk-free interest rates and is contained in options_data.mat
% ***** MAIN INSTRUCTION ******


%% A. SIMULATION TYPE SELECTION
%   Simply comment out the unwanted selection

%  Select a Model
Data.Model_choice = 'md_RBF';   % Radial Basis Function
Data.Model_choice = 'BS';       % Black-Scholes Model

%  Select a filter
Data.Filter_choice = 'SOEKF';   % Second Order Extended Kalman Filter
% Data.Filter_choice = 'EKF';   % Extended Kalman Filter

%  Select type of data to simulate on
Data.Data_type = 'real';        % Real data
% Data.Data_type = 'synthetic';   % Synthetic data (Constant States)

% Select Options contract type
Data.option_type = 'call';      % Call Options contract
% Data.option_type = 'put';     % Put Options contract

%  Select a year for the data and the Strike Price
%  The years available and their respective available strike prices are:
%   2014 - 6800 or 6700; 2013 - 6100 or 6000; 2012 - 5500 or 5600

Data.data_yr = 2014;    % Year
Data.X = 6700;          % Strike Price


%%  B. PARTAMETER INITIALIZATION

if strcmp(Data.Model_choice, 'BS')
    Data.X_a = [(0.1), (0.1)];             % Initial Prior state X[1|1] for synthetic data
    if strcmp(Data.Data_type,'real')
        Data.X_a = [log(0.1), log(0.1)];   % Initial Prior state X[1|1] for real data
    end
    Data.Nx = size(Data.X_a,2);
    Data.P_a = 1e-1*eye(Data.Nx);          % Initial prior state-error covariance P[1|1]
    
    % Set the noise covariances
    Data.Q = 1e-10*eye(Data.Nx);           % State noise covariance Q
    Data.R = 1e-10;                        % Measurement noise covariance R
    
elseif strcmp(Data.Model_choice, 'md_RBF')
    Data.X_a = 1e0*ones(4,2);                    % Initial Prior state X[1|1] for real data
    if strcmp(Data.Data_type,'synthetic')
        Data.X_a = [0.52 0.25]; % 1e0*ones(1,2); % % Initial Prior state X[1|1] for real data
    end
    Data.Nx = size(Data.X_a,1);
    Data.P_a = 1e0*eye(Data.Nx);                % Initial prior state-error covariance P[1|1]
    Data.Q = 1e-6*eye(Data.Nx);                 % State noise covariance Q
    Data.R = 1e-6;                              % Measurement noise covariance R
    
    Data.get_real_means = 0;        % 1- to get the real mean vectors; 0 - otherwise
end

%%  C.  STEADY STATE PARAMETER INITIALIZATION

%  Set the Steady-state Parameters for Synthetic data
if strcmp(Data.Data_type, 'synthetic')
    Data.r = 0.005;      % Constant interest rate of 0.5 percent for BS synthetic data
    Data.sig = 0.1;      % Constant volatility for BS synthetic data @ 10 percent mean deviation
    Data.mn = [0.2 0.5]; % Intialization of constant mean values for RBF synthetic data
end

%%  D.  MODEL AND DATA PREPARATION
%   Generate the model: Output stores the equations for the either
%   Black-Scholes or the RBF depending on the model choice
Model = get_model(Data.Nx, Data.Model_choice);


%  Generate data: Outputs either real or synthetic data that suites the choice
%  of the model
Data = get_data(Model, Data);

%Data.Z_X = Data.C_X;

switch Data.Model_choice
    case 'BS'
        Data.model_in = Data.S;
        Data.model_out = Data.Call;
        if strcmp(Data.option_type,'put')
            Data.model_out = Data.Put;
        end
        Data.states = [Data.sig,Data.r];
        if length(Data.X_a) == 1
            Data.states = [Data.sig];
        end
    case 'md_RBF'
        seq = 1:size(Data.X_a,1);
        Data.w.non_lin_wts = seq./(sum(seq));
        Data.RBFwts = Data.w;
        if strcmp(Data.Data_type, 'synthetic')
            Data.model_in = Data.synthetic_in;
            Data.model_out = Data.synthetic_out;
        elseif strcmp(Data.Data_type, 'real')
            Data.model_in = Data.real_in;
            Data.model_out = Data.C_X;
        end
        Data.states = Data.mean;
end

%%  E.  KALMAN FILTER ITERATION PROCESS

for k = 3: 210 % length(Data.model_in)-3
    t_m = Data.real_in(k,2);        % Time to Maturity
    
    [Data] = KF(Model,Data, k);     % Kalman Filter function
    
    switch Data.Model_choice
        
        
        case 'BS'
            % Assigning posteriori states to new variables
            sig = (Data.X_a(1));
            if strcmp(Data.Data_type,'real')
                sig = exp(Data.X_a(1));
            end
            if Data.Nx == 1
                r = Data.r(k);
            elseif Data.Nx == 2
                r = (Data.X_a(2));
                if strcmp(Data.Data_type,'real')
                    r = exp(Data.X_a(2));
                end
            end
            
            % Options price Estimation
            if strcmp(Data.option_type,'call')
                Price_pred(k,:) = Model.C(Data.S(k),Data.X,r,sig,t_m);
            elseif strcmp(Data.option_type,'put')
                Price_pred(k,:) = Model.P(Data.S(k),Data.X,r,sig,t_m);
            end
            
        case 'md_RBF'
            if strcmp(Data.Data_type,'synthetic')
                Price_pred(k,:) = Model.Mahal_dist(Data.model_in(k,:), Data.X_a, Data.cov);
            elseif strcmp(Data.Data_type,'real')
                Price_pred(k,:) = Model.RBF (Data.model_in(k,:), Data.X_a, Data.cov, Data.RBFwts);
            end
    end
    
    state_est(k,:) = reshape(Data.X_a,[1 size(Data.X_a,1)*size(Data.X_a,2)]);  % Store the states
    eigs(k,:) = max(eig(Data.P_a)');                    % Store the uncertainty measure
    P_mat{k,:} = Data.P_a;                              % Store the state error covariance
    
    sprintf('Busy... Iteration %d of %d ',k,length(Data.model_in))
end

%%  F. PERFORMANCE ANALYSIS

abs_err = abs(Price_pred(1:k,1) - Data.model_out(1:k,1));  % Absolute Error

%  Root Mean Squared Error
switch Data.Model_choice
    case 'BS'
        RMSE = sqrt((mean(Price_pred(101:k,1)./Data.X - Data.model_out(101:k,1)./Data.X).^2))
        RMSE2 = sqrt((mean(Price_pred(5:100,1)./Data.X - Data.model_out(5:100,1)./Data.X).^2))
    case 'md_RBF'
        RMSE = sqrt((mean(Price_pred(101:k,1) - Data.model_out(101:k,1)).^2))
        RMSE2 = sqrt((mean(Price_pred(5:100,1) - Data.model_out(5:100,1)).^2))
end
info = sprintf('Model - %s| Filter - %s| Data - %s',Data.Model_choice,Data.Filter_choice,Data.Data_type)

%Make a plot of the Observations
createPlot([Price_pred(3:k,1),Data.model_out(3:k,1)], 'observation', info, Data)

%make a plot of the states
createPlot([(state_est(3:k,:)), (Data.states(3:k,:))], 'states', info, Data)

hold off;
