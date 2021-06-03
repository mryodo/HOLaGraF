clear all;
close all;
format long;

fileID=1;

edges=readEdges('8.edges');
trigs=readTrigs('8.trigs');

n=8;
G.B1=B1fromEdges(n, edges);
G.B2=B2fromTrig(n, edges, trigs);
G.w=[1 1 1 1 1 1 1 1 1 1]'; G.eps0=0; G.e=0;
%G.w=[1.55 0.86 0.22 1.62 0.3 1.02 0.69 1.87 0.71 1.57]';
G.eps0=0.1;

rng("shuffle");
G.e=-rand(size(G.w));
G.e=G.e/norm(G.e, 2);

k=1;
p=0;
thrs.mu=1.0; thrs.alph=1.; 
thr0=1e-12;

grad=getGrad(G, k, p, thrs, thr0, 1, fileID);
numgrad=getGradNum(G, k, p, thrs, thr0, 1, fileID);

grad-numgrad
