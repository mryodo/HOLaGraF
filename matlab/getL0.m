function L0 = getL0(B1, w, e, eps)
    W=diag(sqrt(w)+eps*e);
    D12=diag(1./(sum(abs(B1)*W, 2)));
    L0=D12*B1*W*W*B1'*D12;
end