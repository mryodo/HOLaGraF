function grad_num=getGradNum(G, k, p, thrs, thr0, normcor, fileID)
    %returns numerical gradient at the point e
    
    myInit=getFk_l2(G, k, p, thrs, fileID);
    
    grad_num=zeros(size(G.e));
    delt=1e-4;
    for i=1:size(G.e, 1)
        if abs(sqrt(G.w(i))+G.eps0*G.e(i))>thr0
            mov=zeros(size(G.e));
            mov(i)=delt;
            G1=G;
            G1.e=G1.e+mov;
            newVal=getFk_l2(G1, k, p, thrs, fileID);
            grad_num(i)=(newVal-myInit)/delt;
        end
    end
    grad_num=grad_num/G.eps0;
    
    if normcor
        matmask = not(abs(sqrt(G.w)+G.e*G.eps0) < thr0);
        PE = G.e.*matmask;
        kappa = dot(grad_num.* matmask, PE)/dot(PE, PE);
        grad_num = grad_num-kappa*PE;
    end

end