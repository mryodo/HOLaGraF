function [eps_fin, e_fin, ix]=getAnswerFromBack(eps_back, es_ans_back, func_fin_back)
    ix=find(func_fin_back>2*1e-5,1)-1;
    eps_fin=eps_back(ix);
    e_fin=es_ans_back(:, ix);
end