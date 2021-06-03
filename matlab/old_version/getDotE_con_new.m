function dE = getDotE_con_new(B1, B2, w, e, eps, k, thr, alph, mu)
    matmask = diag(not(abs(sqrt(w)+e*eps) < thr));
    E = diag(e);
    PE = E.*matmask;

    L1_E = HodgeLW_fr(B1, B2, w, e, eps);
    [U, d]=eig(L1_E,'vector');
    [~, ind] = sort(d);
    U = U(:, ind);
    
    GE = zeros(size(L1_E));
    for i=1:(k+1)
        Gi=getG_i(B1, B2, L1_E, U(:,i), w, e, eps);
        GE = GE + Gi;
    end
    
    L0=getL0(B1, w, e, eps);
    W=diag(sqrt(w)+eps*e);
    D12=diag(1./sqrt(sum(abs(B1)*W, 2)));
    [U, d]=eig(L0,'vector');
    [d, ind] = sort(d);
    U = U(:, ind);
    x = U(:, 2);
    if d(2)<mu
      Add = 2*B1.T*D12*(x*x.T)*D12*B1*W-...
          diag(diag(Sym(B1*W*W*B1.T*D12*(x*x.T))*(D12^3)).T*abs(B1));
      GE = GE - alph/mu*max(0, 1-d(2)/mu)*Add;
    end
    
    kappa = scal(-GE, PE)/scal(PE, PE);
    PGE = GE.*matmask;
    dE = -PGE-kappa*PE;
end