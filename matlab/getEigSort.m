function [d, U]=getEigSort(A)
    [U, d] = eig(A, 'vector');
    [d, ind] = sort(d);
    U = U(:, ind);
end