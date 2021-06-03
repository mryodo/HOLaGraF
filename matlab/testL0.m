clear all;
close all;
format long;

fileID=1;

edges=readEdges('8.edges');
trigs=readTrigs('8.trigs');


n=8;
G.B1=B1fromEdges(n, edges);
G.B2=B2fromTrig(n, edges, trigs);
G.w=[0.2 0.2 1 1 1 1 1 1 1 1]'; G.eps0=0; G.e=0;

figure
tiledlayout(1,4, 'Padding', 'compact', 'TileSpacing', 'compact'); 
nexttile
simpleDrawB1(G.B1, G.w, 0, 0, 0);

L0=getL0(G);
d=eig(L0, 'vector'); d=sort(d); d(2)

G.w=[0.05 0.05 1 1 1 1 1 1 1 1]'; G.eps0=0; G.e=0;
nexttile
simpleDrawB1(G.B1, G.w, 0, 0, 0);

L0=getL0(G);
d=eig(L0, 'vector'); d=sort(d); d(2)

G.w=[0.005 0.005 1 1 1 1 1 1 1 1]'; G.eps0=0; G.e=0;
nexttile
simpleDrawB1(G.B1, G.w, 0, 0, 0);

L0=getL0(G);
d=eig(L0, 'vector'); d=sort(d); d(2)

G.w=[0.0005 0.0005 1 1 1 1 1 1 1 1]'; G.eps0=0; G.e=0;
nexttile
simpleDrawB1(G.B1, G.w, 0, 0, 0);

L0=getL0(G);
d=eig(L0, 'vector'); d=sort(d); d(2)

