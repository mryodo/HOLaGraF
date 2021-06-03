function B2 = B2fromTrig(n, edges, trigs)
    m = size(edges, 1);
    del = size(trigs, 1);
    B2 = zeros(m, del);
    
    for i=1:del
        B2(find(ismember(edges, [trigs(i, 1) trigs(i, 2)], 'rows'), 1), i)=1;
        B2(find(ismember(edges, [trigs(i, 1) trigs(i, 3)], 'rows'), 1), i)=-1;
        B2(find(ismember(edges, [trigs(i, 2) trigs(i, 3)], 'rows'), 1), i)=1;
    end 

end