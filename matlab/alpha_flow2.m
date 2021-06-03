function [e_log, func_log, endings, G, thrs]=alpha_flow2(G, k, alstart, alfinish, thrs, initial, fileID)

    my_beta=1.2; h0=1e-1; thr0=1e-12; 
    
    if ~initial
        G1=G; G1.e=ones(size(G.w))/sqrt(size(G.w, 1)); G1.eps0=1e-4;
        L1=HodgeLW_fr(G1); p=sum(abs(L1) * ones(size(L1, 1), 1)<1e-3);
        thrs.alph=alstart;
        e0=getGrad(G1, k, p, thrs, thr0, 1, fileID); e0=-e0/norm(e0, 2);
        %e0=getGradNumF1(G1, k, p, thrs, thr0, 1, fileID); e0=-e0/norm(e0, 2);
        G.e=e0;
        %rng("shuffle");
        %G.e=-rand(size(G.e));
        %G.e=-[0 0 0 0 0 0 0 0 1 0]'; 
        %G.e=G.e/norm(G.e, 2);
        
    else
        G1=G;
        L1=HodgeLW_fr(G1); p=sum(abs(L1) * ones(size(L1, 1), 1)<1e-3);
        %thrs.alph
        %thrs.alph=alstart;
    end
    %if getFk2(G, k, p, thrs, fileID)>1e-6
        %fprintf("Penalization from the start! Abort, abort!\n");
        %return         
    %end
    
    als=(linspace(sqrt(alstart), sqrt(alfinish), 15)).^2;
    [e_ans, track, ~, ending, ~, ~, ~, ~, G]=single_run(G, k, thrs, h0, my_beta, thr0, fileID);
    %fprintf("alpha: %f  ||  value: %f || ended with: %d \n", thrs.alph, track(end), ending);
    e_log=[e_ans]; func_log=[track(end)]; endings=[ending];
    i_start=2;
    if initial
        i_start=length(als)+1;
    end

    for i=i_start:length(als)
        G.e=e_ans; thrs.alph=als(i);
        [e_ans, track, ~, ending, ~, ~, ~, ~, G]=single_run(G, k, thrs, h0, my_beta, thr0, fileID);
        %fprintf("alpha: %f  ||  value: %f || ended with: %d \n", thrs.alph, track(end), ending);
        e_log=[e_log e_ans]; func_log=[func_log track(end)]; endings=[endings ending];
    end
end