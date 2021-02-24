clear all;
close all;



n=8;
edges=readEdges('8.edges');
%w0=ones(length(edges) , 1);
w0=getRandomWeights(edges);
for i=9:9
    edges=readEdges('8.edges');
    trigs=readTrigs('8.trigs');
    w=w0;
    B1=B1fromEdges(n, edges);
    B2=B2fromTrig(n, edges, trigs);
    %w=ones(length(B1) , 1);
    

    L1=HodgeLW_fr(B1, B2, w, 0, 0);
    L0=getL0(B1, w, 0, 0);
    %L1=B1'*B1+B2*B2';
    %L0=B1*B1';
    spL0=sort(eig(L0));
    spL1=sort(eig(L1));
    spBB=sort(eig(B1'*B1));
    M = zeros(3,max([length(spL0),length(spL1), length(spBB)]));
    M(1,1:length(spL0)) = spL0;
    M(2,1:length(spL1)) = spL1;
    M(3,1:length(spBB)) = spBB;
    
    purge=edges(i, :);
   
    pt=[];
    for j=1:size(trigs, 1)
        if ( (purge(1)==trigs(j, 1)) && ...
                ((purge(2)==trigs(j, 2)) || (purge(2)==trigs(j, 3)) ) ) || ...
                ( (purge(1)==trigs(j, 2)) && ...
                ((purge(2)==trigs(j, 3)) ) ) 
            pt=[pt j];
        end
    end
    trigs(pt, :)=[];
    edges(i, :)=[];
    w(i)=[];
    
    B1=B1fromEdges(n, edges);
    B2=B2fromTrig(n, edges, trigs);
    %w=ones(length(B1) , 1);
    
    thr=1e-8;
    L1=HodgeLW_fr(B1, B2, w, 0, 0);
    L0=getL0(B1, w, 0, 0);
    %L1=B1'*B1+B2*B2';
    %L0=B1*B1';
    
    f=figure('Renderer', 'painters', 'Position', [10 10 450 750], 'visible', 'off');
    subplot(3,1,1)

    W=diag(sqrt(w));
    Aw=getAdjWB1(B1, W);
    G=graph(Aw);
    colors=["#7A76C2"; "#ff6e9c98"; "#f62196"; "#18c0c4"; "#f3907e"; "#66E9EC"];

    plot(G, 'EdgeLabel', sqrt(G.Edges.Weight), 'EdgeColor', colors(4), ...
        'LineWidth', 8*G.Edges.Weight, ...
        'MarkerSize', 20, 'NodeColor', colors(5), ...
        'EdgeFontSize', 16, 'NodeFontSize', 16, 'NodeFontWeight', 'bold')
    title('Graph with eps=', 0, 'fontsize', 20);
    subplot(3, 1, 2)
    
    M_old=M;
    spL0=sort(eig(L0));
    spL1=sort(eig(L1));
    spBB=sort(eig(B1'*B1));
    M = zeros(3,max([length(spL0),length(spL1), length(spBB)]));
    M(1,1:length(spL0)) = spL0;
    M(2,1:length(spL1)) = spL1;
    M(3,1:length(spBB)) = spBB;
    spL0(find(spL0 > thr, 1))-spL1(find(spL1 > thr, 1))
    
    bar(M');
    subplot(3,1,3)
    bar(-(M_old(:, 1:end-1))'+M');
    %print -djpeg test.jpg
    saveas(f, sprintf("test%d.jpg", i));
    close(f);
end


