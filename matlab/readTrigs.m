function trigs = readTrigs(filename)
    opts = detectImportOptions(filename, 'FileType', 'text');
    pre=readmatrix(filename, opts);
    pre=sort(pre, 2);
    [~,idx] = sort(pre(:,1));
    trigs = pre(idx,:);
end

