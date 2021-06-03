function trigs = readTrigs(filename)
    opts = detectImportOptions(filename, 'FileType', 'text');
    pre=readmatrix(filename, opts);
    %size(pre)
    pre=sort(pre, 2);
    %[~,idx] = sort(pre(:,1));
    %trigs = pre(idx,:);
    trigs=sortrows(pre, [1 2 3]);
end

