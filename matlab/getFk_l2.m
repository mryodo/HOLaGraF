function res=getFk_l2(G, k, p, thrs, fileID)
    
    res=getFk1(G, k, p, thrs, fileID)+getFk2(G, k, p, thrs, fileID);
    %+getFk3(G, k, p, thrs, fileID);

end