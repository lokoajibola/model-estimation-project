function bs = Black_Scholes_Model(Nx)
    %BS equations
bs.ln = @(a)log(a)/log(exp(1)); %natural log
bs.N = @(b) normcdf(b);
bs.Np = @(b) normpdf(b);
bs.d_1 = @(S,X,r,sig,t_m)(bs.ln(S/X) + (r + ((sig^2)/2))*sqrt(t_m))/(sig * sqrt(t_m));
bs.d_2 = @(S,X,r,sig,t_m)(bs.ln(S/X) + (r + ((sig^2)/2))*sqrt(t_m))/(sig * sqrt(t_m)) - (sig * sqrt(t_m));
bs.C = @(S,X,r,sig,t_m)(S*bs.N((bs.d_1(S,X,r,sig,t_m)))) - (X*exp(-r*t_m)*bs.N(bs.d_2(S,X,r,sig,t_m)));
bs.P = @(S,X,r,sig,t_m)(S*-bs.N((-bs.d_1(S,X,r,sig,t_m)))) + (X*exp(-r*t_m)*bs.N(-bs.d_2(S,X,r,sig,t_m)));

%First derivatives
bs.dCs = @(S,X,r,sig,t_m) S*sqrt(t_m)*bs.Np(bs.d_1(S,X,r,sig,t_m));
bs.dPs = bs.dCs;
bs.dCr = @(S,X,r,sig,t_m)X * t_m *exp(-r*t_m)*bs.N(bs.d_2(S,X,r,sig,t_m));
bs.dPr = @(S,X,r,sig,t_m)-X * t_m *exp(-r*t_m)*bs.N(-bs.d_2(S,X,r,sig,t_m));

%Second derivatives
%[d^2C/d sigma^2] and [d^2P/d sigma^2]
bs.d2Cs = @(S,X,r,sig,t_m) (S*sqrt(t_m)*bs.d_1(S,X,r,sig,t_m)*bs.d_2(S,X,r,sig,t_m)/sig) *bs.Np(bs.d_1(S,X,r,sig,t_m));
bs.d2Ps = bs.d2Cs;

%[d^2C/d r^2] and [d^2P/d r^2]
bs.d2Cr = @(S,X,r,sig,t_m) -X*t_m*exp(-r*t_m)*(t_m*bs.N(bs.d_2(S,X,r,sig,t_m)) - ((bs.d_2(S,X,r,sig,t_m)*sqrt(t_m)/sig)*bs.Np(bs.d_2(S,X,r,sig,t_m))));
bs.d2Pr = @(S,X,r,sig,t_m) X*t_m*exp(-r*t_m)*(t_m*bs.N(-bs.d_2(S,X,r,sig,t_m)) - ((bs.d_2(S,X,r,sig,t_m)*sqrt(t_m)/sig)*bs.Np(-bs.d_2(S,X,r,sig,t_m))));

%[d^2C/d r d sigma] and [d^2P/d r d sigma]
bs.d2Csr = @(S,X,r,sig,t_m)-S*bs.d_1(S,X,r,sig,t_m)*t_m/sig*(bs.Np(bs.d_1(S,X,r,sig,t_m)));
bs.d2Psr = bs.d2Csr;

%Jacobian and Hessian
if Nx == 2
    bs.Jac1 = @(S,X,r,sig,t_m)[bs.dCs(S,X,r,sig,t_m), bs.dCr(S,X,r,sig,t_m)];
    bs.Jac3 = @(S,X,r,sig,t_m)[bs.dPs(S,X,r,sig,t_m), bs.dPr(S,X,r,sig,t_m)];
    bs.Jac2 = @(S,X,r,sig,t_m)[bs.dCs(S,X,r,sig,t_m), bs.dCr(S,X,r,sig,t_m); bs.dPs(S,X,r,sig,t_m), bs.dPr(S,X,r,sig,t_m)];
    bs.Hes1 = @(S,X,r,sig,t_m)[bs.d2Cs(S,X,r,sig,t_m), bs.d2Csr(S,X,r,sig,t_m); bs.d2Csr(S,X,r,sig,t_m), bs.d2Cr(S,X,r,sig,t_m)];
    bs.Hes2 = @(S,X,r,sig,t_m)[bs.d2Ps(S,X,r,sig,t_m), bs.d2Psr(S,X,r,sig,t_m); bs.d2Psr(S,X,r,sig,t_m), bs.d2Pr(S,X,r,sig,t_m)];
elseif Nx == 1
    bs.Jac1 = @(S,X,r,sig,t_m) bs.dCs(S,X,r,sig,t_m);
    bs.Jac3 = @(S,X,r,sig,t_m) bs.dPs(S,X,r,sig,t_m);
    bs.Hes1 = @(S,X,r,sig,t_m) bs.d2Cs(S,X,r,sig,t_m);
    bs.Hes2 = @(S,X,r,sig,t_m) bs.d2Ps(S,X,r,sig,t_m);
end
bs.choice = 'BS';
end