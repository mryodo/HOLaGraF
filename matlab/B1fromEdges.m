function B1 = B1fromEdges(n, edges)
    m = size(edges, 1);
    B1 = zeros(n, m);
    
    for i=1:m
       B1(edges(i, 1), i) = -1;
       B1(edges(i, 2), i) = 1;
    end
end