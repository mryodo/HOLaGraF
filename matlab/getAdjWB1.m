function A = getAdjWB1(B1, W)   
    A = diag(diag(B1*W*W*B1'))-B1*W*W*B1';
end