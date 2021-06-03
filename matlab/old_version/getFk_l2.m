function res=getFk_l2(L1, k)
    d=eig(L1,'vector');
    d = sort(d);
    res=sqrt(sum(power(d(1:k+1), 2)));
end