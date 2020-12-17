function simpleDrawB1(B1, w, eps, e)
    W=diag(sqrt(w)+eps*e);
    Aw=getAdjWB1(B1, W);
    G=graph(Aw);
    colors=["#7A76C2"; "#ff6e9c98"; "#f62196"; "#18c0c4"; "#f3907e"; "#66E9EC"];
    figure
   
    plot(G, 'EdgeLabel', sqrt(G.Edges.Weight), 'EdgeColor', colors(4), ...
        'LineWidth', 8*G.Edges.Weight, ...
        'MarkerSize', 20, 'NodeColor', colors(5), ...
        'EdgeFontSize', 16, 'NodeFontSize', 16, 'NodeFontWeight', 'bold')
     title('Graph with eps=', eps, 'fontsize', 20);

end