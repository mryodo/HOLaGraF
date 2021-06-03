clear all;
close all;



n=8;
edges=readEdges('8.edges');
trigs=readTrigs('8.trigs');

B1=B1fromEdges(n, edges);
B2=B2fromTrig(n, edges, trigs);

