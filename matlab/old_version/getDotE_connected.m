function dE = getDotE_connected(B1, B2, w, e, eps, k, thr, alpha)
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
    [U, d]=eig(L0,'vector');
    [~, ind] = sort(d);
    U = U(:, ind);
    x = U(:, 2);
    GE = GE + 2*alpha*B1'*(x*x')*B1*diag(sqrt(w)+eps*e);
    
    kappa = scal(-GE, PE)/scal(PE, PE);
    PGE = GE.*matmask;
    dE = -PGE-kappa*PE;
end