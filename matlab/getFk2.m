function res=getFk2(G, k, p, thrs, fileID)
    alph=thrs.alph;
    mu=thrs.mu;
    
    L0=getL0(G);

    d=eig(L0,'vector');
    d=sort(d);
    res=0.5*alph*(max(0, 1-d(2)/mu))^2;
end