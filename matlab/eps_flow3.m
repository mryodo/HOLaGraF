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
G.w=[1.55 0.86 0.22 1.62 0.3 1.02 0.69 1.87 0.71 1.57]';

k=1;
L0=getL0(G); tmp=sort(eig(L0));

thrs.mu=tmp(2)*1.0; thrs.alph=0.; 
alstart=1; alfinish=100;

%eps_s=linspace(0.9*min(sqrt(G.w)), 2.5, 20);
Heps=0.025;
Heps2=Heps/2;

func_fin_forw=[]; es_ans_forw=[]; eps_forw=[];
func_fin_back=[]; es_ans_back=[]; eps_back=[];

while 1
    if length(func_fin_forw)>0
        if func_fin_forw(end)<1e-5
            break
        end
    end
    if G.eps0>2.5
        break
    end
    G.eps0=G.eps0+Heps;
    if G.eps0==Heps
        [e_log, func_log, endings]=alpha_flow2(G, k, alstart, alfinish, thrs, 0, fileID);
    else
        [e_log, func_log, endings]=alpha_flow2(G, k, alstart, alfinish, thrs, 1, fileID);
    end
    fprintf("eps: %f   |||   func: %f  \n", G.eps0, func_log(end));
    G.e=e_log(:, end);
    func_fin_forw=[func_fin_forw func_log(end)];
    es_ans_forw=[es_ans_forw e_log(:, end)];
    eps_forw=[eps_forw G.eps0];
end

while 1
    if length(func_fin_back)>0
        if func_fin_back(end)>1e-2
            break
        end
    end
    if G.eps0<=0.05
        break
    end
    G.eps0=G.eps0-Heps2;
    [e_log, func_log, endings]=alpha_flow2(G, k, alstart, alfinish, thrs, 1, fileID);
    
    fprintf("eps: %f   |||   func: %f  \n", G.eps0, func_log(end));
    G.e=e_log(:, end);
    func_fin_back=[func_fin_back func_log(end)];
    es_ans_back=[es_ans_back e_log(:, end)];
    eps_back=[eps_back G.eps0];
end


figure
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact'); 
nexttile
hold on
plot(eps_forw, func_fin_forw, 'linewidth', 4);
plot(eps_back, func_fin_back, 'linewidth', 4);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Epsilons', 'fontsize', 16);
ylabel('Final Functional', 'fontsize', 16);
nexttile
simpleDrawB1(G.B1, G.w, 0, 0, 0);
nexttile
simpleDrawB1(G.B1, G.w, eps_forw(end), es_ans_forw(: ,end), 0);

nexttile
%lab={'[1, 2]', '[1, 3]', '[1, 4]', '[2, 3]', '[3, 4]'};
h1=heatmap([es_ans_forw es_ans_back]);
%h1.YDisplayLabels=lab;
