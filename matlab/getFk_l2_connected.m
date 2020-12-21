function res=getFk_l2_connected(L1, k, L0, alpha)
    d=eig(L1,'vector');
    d = sort(d);
    res=sum(power(d(1:k+1), 2));
    d=eig(L0,'vector');
    res=res+alpha*d(2)^2;
end