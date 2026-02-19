function bs = BS_Model1()

syms  S X r sig t_m

bs.N = @(b) normcdf(b);
bs.Np = @(b) normpdf(b);
bs.d_1 = @(S,X,r,sig,t_m)(log(S/X) + (r + ((sig^2)/2))*sqrt(t_m))/(sig * sqrt(t_m));
bs.d_2 = @(S,X,r,sig,t_m)(log(S/X) + (r + ((sig^2)/2))*sqrt(t_m))/(sig * sqrt(t_m)) - (sig * sqrt(t_m));

%Black-Scholes Call and Put options equation
bs.C = (S*bs.N((bs.d_1(S,X,r,sig,t_m)))) - (X*exp(-r*t_m)*bs.N(bs.d_2(S,X,r,sig,t_m)));
bs.P = (-S*bs.N((-bs.d_1(S,X,r,sig,t_m)))) - (X*exp(-r*t_m)*bs.N(-bs.d_2(S,X,r,sig,t_m)));

bs.choice = 'BS';
end