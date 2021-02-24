function dE = getDotE_con_new_diag(B1, B2, w, e, eps, k, thr, alph, mu)
    matmask = not(abs(sqrt(w)+e*eps) < thr);
    PE = e.*matmask;

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
    D12=diag(1./sum(abs(B1)*W, 2));
    [U, d]=eig(L0,'vector');
    [d, ind] = sort(d);
    U = U(:, ind);
    x = U(:, 2);
    if d(2)<mu
      Add = 2*B1'*D12*(x*x')*D12*B1*W-...
          2*diag(diag(Sym(D12*(x*x')*L0))'*abs(B1));
      GE = GE - alph/mu*max(0, 1-d(2)/mu)*Add;
    end
    GE=diag(GE);
    
    kappa = dot(GE, PE)/dot(PE, PE);
    PGE = GE.*matmask;
    dE = -PGE+kappa*PE;
end