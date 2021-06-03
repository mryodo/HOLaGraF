function res=getFk3(G, k, p, thrs, fileID)
    C2=thrs.C2;

    L1=HodgeLW_fr(G);

    [U, d]=eig(L1, 'vector');
    [d, ind] = sort(d);
    U = U(:, ind);

    res=0;
    for ii=1:size(L1, 1)
       x=U(:, ii);
       for j=1:size(L1, 1)
           if ( sqrt(G.w(j))+G.eps0*G.e(j) ) <= 0.1*sqrt(G.w(j))
               res=res + 0.5*C2*(max( abs(x(j))/0.8 - 1, 0))^2;
               %fprintf(fileID, "versor at %d \n", j);
           end
       end
    end
end