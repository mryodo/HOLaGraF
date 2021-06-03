function res=getFk_l2_con_new(L1, k, L0, alph, mu)
%%% We can improve efficiency here by computing only the eigevects we need
    d=eig(L1,'vector');
    d = sort(d);
    res=sum(power(d(1:k+1), 2));
    d=eig(L0,'vector');
    d=sort(d);
    res=0.5*res+0.5*alph*(max(0, 1-d(2)/mu))^2;
end