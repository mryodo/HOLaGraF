function M=getM(G, W, Dt)
    M=zeros(size(G.B2'));
    for i=1:size(M, 1)
        ix=find( abs(diag(W)-Dt(i,i))<1e-08, 1);
        M(i, ix)=1;
    end
end