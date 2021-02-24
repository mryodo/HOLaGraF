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
mu=tmp(2)*0.9;
%w=[1.01 1.02 .99 0.98 1.08 1.05 0.95 0.25 1.1 0.99]';

simpleDrawB1(B1, w, 0, 0);


e=-ones(size(w));
e=e/norm(e, 2);
alph=10.;

k=1;
h=1e-1;
beta=1.2;
k=1;
eps=sqrt(w(8));


eps=0;
e=ones(size(w));
e=e/norm(e, 2);
e=getDotE_con_new_diag(B1, B2, w, e, eps, k, thr, alph, mu);
e=-e/norm(e, 2);

eps=sqrt(w(8));

display("initial setup done");

L1_E = HodgeLW_fr(B1, B2, w, e, eps);
L0=getL0(B1, w, e, eps);
track = [getFk_l2_con_new(L1_E, k, L0, alph, mu)];

t_cur = 0;
log = [];
ts = [0];
h0=h;
tfin=1000;
nu=1.0;
i=0;

minval=1.0;
emin=zeros(size(e));

h_log=[];
h_desc=[];

for i=1:10000 
    e0 = e;
    dE = getDotE_con_new_diag(B1, B2, w, e, eps, k, thr, alph, mu);
    while 1
        e1 = e-h*dE;
        e1(sqrt(w)+eps*e1<0)=-1./eps*sqrt(w(sqrt(w)+eps*e1<0));
        e1=e1/norm(e1, 2);
        L1_E = HodgeLW_fr(B1, B2, w, e1, eps);
        L0=getL0(B1, w, e1, eps);
        newval=getFk_l2_con_new(L1_E, k, L0, alph, mu);
        if  (newval < track(end))      
            fprintf("%d / step %f accepted\n", i, h);
            h=h*beta;
            h_log=[h_log h];
            h_desc=[h_desc 1];
            e = e1;
            break
        else
            e=e0;
            fprintf("%d / step %f rejected\n", i, h);
            h=h/beta;
            h_log=[h_log h];
            h_desc=[h_desc 2];
        end
        if h<1e-10
            e=e+ nu*getDotE_con_new_diag(B1, B2, w, e, 0, k, thr, alph, mu);
            e=e/norm(e, 2);
            h=h0;
            fprintf("%d / step %f resetted\n", i, h);
            break
        end
        
    end
    log=[log sort(eig(L1_E))'];
    track=[track newval];
    if track(end)<minval
        minval=track(end);
        emin=e;
    end
    t_cur = t_cur + h;
    ts=[ts t_cur];
end



figure
p1=loglog(ts, (track), 'linewidth', 4);
p1.Color(4)=0.5;

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time', 'fontsize', 16);
ylabel('Functional', 'fontsize', 16);

simpleDrawB1(B1, w, eps, e);
simpleDrawB1(B1, w, eps, emin);
