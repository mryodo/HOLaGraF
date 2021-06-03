function L1=HodgeLW_fr(G)
    W = diag(sqrt(G.w)+G.eps0*G.e);
    Dt = getDt(G.B2, W);
    L1 = W*(G.B1'*G.B1)*W+inv(W)*G.B2*Dt*Dt*G.B2'*inv(W);
    
    % It produces symmetric up to 1e-15 matrix. Need to be better?
end