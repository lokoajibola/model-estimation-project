function [Data] = get_data(Model, Data)

%Load in available asset prices, options prices and interest rates
load('options_data.mat')


if strcmp(Data.Data_type, 'synthetic')
    r = Data.r;
    sig = Data.sig;
    mn = Data.mn;
end

%  Assignment of data to variable   
switch Data.Data_type

    case 'real'
        data_yr = Data.data_yr;
        switch data_yr
            case 2014
                if Data.X == 6700
                    loaded_data = optn.d14_6700; 
                elseif Data.X == 6800
                    loaded_data = optn.d14_6800;
                end
                S = loaded_data(:,1);
                Data.Call = loaded_data(:,2);
                Data.Put = loaded_data(:,3);
                rate = optn.rates_14(1:length(S),3);
            case 2013
                if Data.X == 6000
                    loaded_data = optn.d13_6000; 
                elseif Data.X == 6100
                    loaded_data = optn.d13_6100;
                end
                S = loaded_data(:,1);
                Data.Call = loaded_data(:,2);
                Data.Put = loaded_data(:,3);
                rate = optn.rates_13(1:length(S),3);
            case 2012
                if Data.X == 5600
                    loaded_data = optn.d12_5600;  
                elseif Data.X == 5500
                    loaded_data = optn.d12_5500;
                end
                S = loaded_data(:,1);
                Data.Call = loaded_data(:,2);
                Data.Put = loaded_data(:,3);
                rate = optn.rates_12(1:length(S));
        end
        
        Data.S = S;             % Price of underlying asset(FTSE100)
        X = Data.X;             % Strike price
        Data.r = rate./100;     % Interest rates
        Data.T = length(S);     % Total Time to maturity of option
        Data.yr = 253;          % Total number of days in a year as seen in Hutchinson
        Data.C_X = Data.Call./Data.X; % Call option prices scaled down by the strike price
        Data.P_X = Data.Put./Data.X;  % Put option prices scaled down by the strike price
        
    case 'synthetic'
        S = optn.d14_6700(:,1); % Price of underlying asset(FTSE100)
        X = 6450;               % Strike price
        Data.S = S;
        Data.X = X;
        Data.r = repmat(r,length(S),1); % Risk-free interest rates
        Data.sig = repmat(sig,length(S),1); %Volatility
        Data.T = length(S); % Total time to maturity of option
        Data.yr = 253;      % Total days in a year
        
        % Generate the option prices using constant interest and volatility
        for i = 1:Data.T
            t_m = (Data.T - i)/Data.yr;
            [C1(i,:), P1(i,:)] = blsprice(S(i), X, Data.r(i), t_m, Data.sig(i));
            Data.Call(i,:) = C1(i,:) + abs(normrnd(0,0.001));   % Options prices
            Data.C_X(i,:) = C1(i,:)/X + abs(normrnd(0,0.00001));% Scaled options prices
            Data.Put(i,:) = P1(i,:) + abs(normrnd(0,0.001));    % Options prices
            Data.P_X(i,:) = P1(i,:)/X + abs(normrnd(0,0.00001));% Scaled options prices
        end
        
end

%The real input data for the RBF model
Data.real_in = [S./X, (Data.T - (1:length(S))')./Data.yr];

switch Data.Model_choice
    case 'BS'
        % Generate the Historical and implied volatilities
        for i = 1:Data.T
            t_m = (Data.T - i)/Data.yr;
            Data.t_m(i,:) = t_m;
            %Historical volatility
            days = 50; % Window for calculating historical volatility
            if i < days+1
                Data.hist_vol(i,:) = std(price2ret(S(1:days)));
            else
                Data.hist_vol(i,:) = std(price2ret(S(i-days:i)));
            end
            
            %Implied volatility
            Data.imp_vol(i,:) = blsimpv(S(i), X, Data.r(i), t_m, Data.Call(i,:));
            Data.imp_vol2(i,:) = blsimpv(S(i), X, Data.r(i), t_m, Data.Put(i,:));
            
        end
        
        Data.sig = Data.imp_vol;
        if strcmp(Data.option_type,'put')
            Data.sig = Data.imp_vol2;
        end
        
    % Generate data for RBF model
    case 'md_RBF'
        
        %Generate the real mean vectors for the real options data
        if strcmp(Data.Data_type,'real')
            % Initialization of covariance to identity matrix
            Data.cov = eye(2);
            mn = []; new_mn = [];
            Data.w.bias = -0.2;
            Data.w.lin_wts = [-0.1, -0.1];
            if Data.get_real_means == 1
                for i = 2:length(Data.real_in)
                    [mn(i,:), new_mn] = get_means1(Data.real_in(i,:), Data.w, Data.C_X(i,:),size(Data.X_a,1), new_mn,i);
                    i
                end
                Data.mean = mn;
            else
                Data.mean = zeros(300,2);
            end
        elseif strcmp(Data.Data_type,'synthetic')
            
            mn = [0.2 0.5]; % Intialization of constant mean values
            Data.cov = eye(length(mn)); %Initialization of covariance matrix
            
            % Generate 100 random input vectors
            Data.synthetic_in = Data.real_in;
            Data.mean = repmat(mn,length(Data.synthetic_in),1);
            % Generate synthesized output data using the RBF mahalanobis distance
            % model
            for i = 1: size(Data.synthetic_in,1)
                
                Data.synthetic_out(i,:) = Model.Mahal_dist(Data.synthetic_in(i,:),Data.mean(i,:),Data.cov) ...
                    + 0.005*randn;
                
            end
        end
        
        
end