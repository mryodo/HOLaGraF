function L0 = getL0(B1, w, e, eps)
    W=diag(sqrt(w)+eps*e);
    L0=B1*W*W*B1';
end