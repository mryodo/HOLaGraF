function B=myinv(A)
    w=diag(A);
    B=diag(1./w);
end