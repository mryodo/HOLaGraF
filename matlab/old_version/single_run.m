clear all;
close all;

edges=readEdges('6.edges');
trigs=readTrigs('6.trigs');

n=7;
B1=B1fromEdges(n, edges);
B2=B2fromTrig(n, edges, trigs);
w=getRandomWeights(edges);

%simpleDrawB1(B1, w, 0, 0);
thr=1e-8;
L1=HodgeLW_fr(B1, B2, w, 0, 0);

%w=[1.01 1.02 .99 0.98 1.08 1.05 0.95 0.25 1.1 0.99]';

simpleDrawB1(B1, w, 0, 0);

e = -rand(size(w));
e = e/norm(e, 2);
eps = sqrt(w(8));

L1_E = HodgeLW_fr(B1, B2, w, e, eps);
track = [getFk_l2(L1_E, k)];
h = 1e-5;
t_cur = 0;
log = [];
ts = [0, ];
beta=2;

for i=1:100
    e0 = e;
    while 1
        while 1
            e = e0;
            dE = getDotE(B1, B2, w, e, eps, 1, thr);
            E1 = diag(e)+h*dE;
            e = diag(E1);
            fprintf('h:  %s  |  time: %f ||   E_norm: %f  , dE-orth: %f ||   F=%f\n', ...
                    h, t_cur,round(scal(E1, E1), 3), round(scal(dE, E1), 3),track(end));
            if sum(sqrt(w)+eps*e < 0) > 0
                h = h/2;
            else
                break
            end
            if h<1e-10
                h=1e-5;
                break
            end
        end
        L1_E = HodgeLW_fr(B1, B2, w, e, eps);
        if getFk_l2(L1_E, 1) > track(end)
            h = h/beta;
        else
            h = h*beta;
            break
        end
    end
    L1_E = HodgeLW_fr(B1, B2, w, e, eps);
    log=[log, sort(eig(L1_E, 'vector'))];
    track=[track, getFk_l2(L1_E, 1)];
    dE = getDotE(B1, B2, w, e, eps, 1, thr);
    t_cur =t_cur + h;
    ts=[ts t_cur];
end

simpleDrawB1(B1, w, e, eps)

figure
plot(ts, track, 'linewidth', 4)
