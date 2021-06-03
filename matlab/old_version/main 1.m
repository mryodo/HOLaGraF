edges=readEdges('6.edges');
trigs=readTrigs('6.trigs');

n=7;
B1=B1fromEdges(n, edges);
B2=B2fromTrig(n, edges, trigs);
w=getRandomWeights(edges);

%simpleDrawB1(B1, w, 0, 0);

L1=HodgeLW_fr(B1, B2, w, 0, 0);
[U, d]=eig(L1,'vector');
[d, ind] = sort(d);
U = U(:, ind)
d

