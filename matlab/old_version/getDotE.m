function dE = getDotE(B1, B2, w, e, eps, k, thr)
    matmask = diag(not(abs(sqrt(w)+e*eps) < thr));
    E = diag(e);
    PE = E.*matmask;

    L1_E = HodgeLW_fr(B1, B2, w, e, eps);
    [U, d]=eig(L1_E,'vector');
    [d, ind] = sort(d);
    U = U(:, ind);
    
    GE = zeros(size(L1_E));
    for i=1:(k+1)
        Gi=getG_i(B1, B2, L1_E, U(:,i), w, e, eps);
        GE = GE + Gi;
    end
    kappa = scal(-GE, PE)/scal(PE, PE);
    PGE = GE.*matmask;
    dE = -PGE-kappa*PE;
end