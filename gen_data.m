% generate synthetic data for log normal S
t = 253;
S = zeros(t,1);
S(1) = 500; 
z = normrnd(0.005/t, 0.01/t, [t 1]);

for i = 2:253
    S(i) = S(1)*exp(sum(z(1:i)));
    
end


%for lognormal distribution
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
X = lognrnd(mu,sigma,1,250);
