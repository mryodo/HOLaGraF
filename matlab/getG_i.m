function res = getG_i(B1, B2, L1_E, x, w, e, eps)
    W = diag(sqrt(w)+eps*e);
    Dt = getDt(B2, W);
    P2 = myinv(W)*B2*Dt*Dt*B2'*myinv(W);
    w_t = diag(W);
    on = ones(size(B2, 2), 1);
    tmp = (on*w_t').*abs(B2');
    M = zeros(size(B2'));
    tmp(tmp == 0) = 2*max(w_t);
    [~, argmax] = max(-tmp, [], 2);
    M(1:size(M,1), argmax) = 1;
    Add = diag(M'*diag(B2'*myinv(W)*x*x'*myinv(W)*B2*Dt));
    Gi = 2.*Sym(B1'*B1*W*x*x'-myinv(W)*x*x'*P2)+2.*Add;
    res = (x'*L1_E*x)*Gi;
end
