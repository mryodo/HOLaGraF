function Dt = getDt(B2, W)
    w=diag(W);
    on=ones(size(B2, 2), 1);
    tmp=(on*w').*abs(B2');
    minval=zeros(size(tmp, 1), 1);
    for i=1:size(tmp, 1)
       if nnz(tmp(i, :)) >= 3
          minval(i)=min(nonzeros(tmp(i, :))); 
       end
    end
    Dt=diag(minval);
end