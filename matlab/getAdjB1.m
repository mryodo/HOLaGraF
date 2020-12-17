function A = getAdjB1(B1)
     A = diag(diag(B1*B1'))-B1*B1';
end