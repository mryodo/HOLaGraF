function grad=getGrad(G, k, p, thrs, thr0, normcor, fileID)
    L1=HodgeLW_fr(G);
    [vals, vecs]=getEigSort(L1);
    W=diag(sqrt(G.w)+G.eps0*G.e);
    Dt=getDt(G.B2, W);
    P2=inv(W) * G.B2 * Dt * Dt * G.B2' * inv(W);
    preMat=G.B1' * G.B1 * W;
    preMat2=Dt * G.B2' * inv(W);
    preMat3=inv(W) * G.B2;
    M=getM(G, W, Dt);

    grad=zeros(size(G.e, 1), 1);
    for i=1:k+p+1
        x=vecs(:, i);
        lam=vals(i);
        dE1=2*lam*Sym( preMat * (x*x') - P2 * (x*x') * inv(W) );
        A=preMat2 * (x*x') *preMat3;
        dE2=2*lam*diag(M' * diag(A));
        grad = grad+diag(dE1+dE2);
    end

    L0=getL0(G);
    [vals, vecs]=getEigSort(L0);
    lam=vals(2);
    x=vecs(:, 2);

    if (1-lam/thrs.mu)>0
        A = getAdjWB1(G.B1, W);
        A = A+eye(size(A, 1));
        d = A * ones(size(A, 1), 1);
        C=-thrs.alph/thrs.mu*max(0, 1-lam/thrs.mu);
        tmp= ((W*G.B1').*(G.B1')) * diag(Sym(diag(d.^(-1))* (x*x')*L0));
        dE4=-2*diag(tmp);
        dE5=2*G.B1'*diag(d.^(-1/2))*(x*x')*diag(d.^(-1/2))*G.B1*W;

        grad=grad+C*diag(dE4+dE5);
    end


    if normcor
        matmask = not(abs(sqrt(G.w)+G.e*G.eps0) < thr0);
        PE = G.e.*matmask;
        kappa = dot(grad.* matmask, PE)/dot(PE, PE);
        grad = grad-kappa*PE;
    end
end