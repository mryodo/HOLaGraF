function res=getFk1(G, k, p, thrs, fileID)
    L1=HodgeLW_fr(G);

    d=eig(L1,'vector');
    d = sort(d);
    res=0.5*sum(power(d(1:p+k+1), 2));
end