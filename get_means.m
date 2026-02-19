function final_mn = get_means(x, cov, out, mn_sz)

bias = 0.001;
wt1 = [0.01, 0.02];
lamda = repmat(0.01, 1, mn_sz);
residual = out - bias - x*wt1';

for i = 1:mn_sz
    z = sym('z', [1 2]);
    residual2 = residual*(lamda(i)/sum(lamda));
    [a,b] = vpasolve(z*cov*z' == residual2,z);
    try
        mn = x - [double(a),double(b)];
    catch
        mn = NaN(1,2);
%         mn = [0,0];
    end
    new_mn(i,:) = mn;
end

final_mn = reshape(new_mn,1,[]);


end