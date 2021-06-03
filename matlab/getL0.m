function L0 = getL0(G)
    W=diag(sqrt(G.w)+G.eps0*G.e);
    A = getAdjWB1(G.B1, W);
    A = A+eye(size(A, 1));
    d = A * ones(size(A, 1), 1);
    L0=diag(d)-A;
    L0=diag(d.^(-1/2))*L0*diag(d.^(-1/2));
    %D12=sqrt(diag(1./(sum(abs(G.B1)*W*W, 2))));
    %L0=D12*G.B1*W*W*G.B1'*D12;
    %L0=G.B1*W*W*G.B1';
end