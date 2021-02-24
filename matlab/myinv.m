function B=myinv(A)
    w=diag(A);
    thr=1e-12;
    tmp=zeros(size(w));
    tmp(w>thr)=1./(w(w>thr));
    B=diag(tmp);
end