function [out] = eval_jac(a, m, c, func)
% This function simply substitutes the components:
% input vector [a]
% mean vector [m]
% covariance [c]
% into the equation of the Jacobian or Hessian in [func]
% and evaluates the result as [out]

switch length(m)
    case 1
        in1 = a(1,1);
        cov11 = c(1,1);
        mn1 = m(1,1);
        
        syms a11 m1  c1_1
        out = eval(subs(func, [a11,m1,c1_1],[in1,mn1,cov11]));
        
    case 2
        in1 = a(1,1);
        cov11 = c(1,1);
        mn1 = m(1,1);
        in2 = a(1,2);
        mn2 = m(1,2);
        cov12 = c(1,2);
        cov21 = c(2,1);
        cov22 = c(2,2);
        
        syms a11 a12 m1 m2 c1_1 c1_2 c2_1 c2_2
        out = eval(subs(func, [a11,a12,m1,m2,c1_1,c1_2,c2_1,c2_2], ...
            [in1,in2,mn1,mn2,cov11,cov12,cov21,cov22]));
        
        
     case 4
        in1 = a(1,1);
        cov11 = c(1,1);
        mn1 = m(1);
        in2 = a(1,2);
        mn2 = m(2);
        mn3 = m(3);
        mn4 = m(4);
        cov12 = c(1,2);
        cov21 = c(2,1);
        cov22 = c(2,2);
        
        syms a11 a12 m1 m2 c1_1 c1_2 c2_1 c2_2 m3 m4
        out = eval(subs(func, [a11,a12,m1,m2,m3,m4,c1_1,c1_2,c2_1,c2_2], ...
            [in1,in2,mn1,mn2,mn3,mn4,cov11,cov12,cov21,cov22]));
end