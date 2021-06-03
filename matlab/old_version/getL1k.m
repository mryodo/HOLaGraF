function res=getL1k(L1, thr)
    res=sum(abs(eig(L1,'vector'))<thr);
end