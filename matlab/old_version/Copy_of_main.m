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
L0=getL0(B1, w, 0, 0);
tmp=sort(eig(L0));
mu=tmp(2)*0.5;
%w=[1.01 1.02 .99 0.98 1.08 1.05 0.95 0.25 1.1 0.99]';

simpleDrawB1(B1, w, 0, 0);

figure
hold on
for rep=1:1
    %eps=1e-8;
    e=-ones(size(w));
    e=e/norm(e, 2);
    alph=10.;
    %mu=0.5;
    k=1;
    %e=diag(getDotE_con_new(B1, B2, w, e, eps, k, thr, alph,mu));
    %e=-e/norm(e,2);
    %e = -rand(size(w));
    h=1e-1;
    beta=1.2;
    k=1;
    eps=sqrt(w(8));
    display("initial setup done");

    %[e, L1_E, ts, track, log]=run(B1, B2, w, eps, e, beta, h, k, thr, 0, 0);
    [e, L1_E, ts, track, log]=run_con_new(B1, B2, w, eps, e, beta, h, k, thr, alph, mu, 0, 0);
    p1=loglog(ts, (track), 'linewidth', 4);
    p1.Color(4)=0.5;
    
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time', 'fontsize', 16);
ylabel('Functional', 'fontsize', 16);

simpleDrawB1(B1, w, eps, e);
