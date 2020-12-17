function L1=HodgeLW_fr(B1, B2, w, e, eps)
    W = diag(sqrt(w)+eps*e);
    Dt = getDt(B2, W);
    L1 = W*B1'*B1*W+myinv(W)*B2*Dt*Dt*B2'*myinv(W);
    
    % It produces symmetric up to 1e-15 matrix. Need to be better?
end