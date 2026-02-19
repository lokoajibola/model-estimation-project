function pred = get_C_X(x, w, mn, cov)
%The weights are contained in w and are preset to constant values with cues
%from hutchinson
for i = 1:size(mn,1)
%     cov1 = cov{i};
    cov1 = eye(2);
    nl(i,:) = sqrt((x - mn(i,:)) * inv(cov1) * (x - mn(i,:))');
end

pred =  w.non_lin_wts*nl + w.lin_wts*x' + w.bias;

end