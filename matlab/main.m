clear all;
close all;

edges=readEdges('6.edges');
trigs=readTrigs('6.trigs');

n=7;
B1=B1fromEdges(n, edges);
B2=B2fromTrig(n, edges, trigs);
w=getRandomWeights(edges);

thr=1e-8;
L1=HodgeLW_fr(B1, B2, w, 0, 0);

%w=[1.01 1.02 .99 0.98 1.08 1.05 0.95 0.25 1.1 0.99]';

simpleDrawB1(B1, w, 0, 0);
min=1.0;
e_min=zeros(size(w));

figure
hold on
for rep=1:7
    e = -rand(size(w));
    h=1e-5;
    beta=2.0;
    k=1;
    eps=sqrt(w(8));

    [e, L1_E, ts, track, log]=run(B1, B2, w, eps, e, beta, h, k, thr, 0, 0);
    p1=loglog(ts, track, 'linewidth', 4);
    p1.Color(4)=0.5;
    
    if track(end)<min
       e_min=e; 
       min=track(end);
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time', 'fontsize', 16);
ylabel('Functional', 'fontsize', 16);

simpleDrawB1(B1, w, eps, e);
