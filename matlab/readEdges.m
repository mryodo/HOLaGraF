function edges=readEdges(filename)
    opts = detectImportOptions(filename, 'FileType', 'text');
    pre=readmatrix(filename, opts);
    pre=sort(pre, 2);
    %[~,idx] = sort([pre(:,1)]);
    %edges = pre(idx,:);
    edges=sortrows(pre, [1 2]);
end