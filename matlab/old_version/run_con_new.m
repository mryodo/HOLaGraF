function [e, L1_E, ts, track, log] = run_con_new(B1, B2, w, eps, e, beta, h, k, thr, alpha, mu, verbose, draw)
    e = e/norm(e, 2);
    L1_E = HodgeLW_fr(B1, B2, w, e, eps);
    L0=getL0(B1, w, e, eps);
    track = [getFk_l2_con_new(L1_E, k, L0, alpha, mu)];
    t_cur = 0;
    log = [];
    ts = [0, ];
    h0=h;
    tfin=10;
    nu=0.1;
    i=1;
    while (ts(end)<tfin) || (i<=100000)
        e0 = e;
        while 1
            while 1
                e = e0;
                dE = getDotE_con_new(B1, B2, w, e, eps, k, thr, alpha, mu);
                E1 = diag(e)-h*dE;
                e = diag(E1);
                e = e/norm(e, 2);
                if verbose
                    fprintf('h:  %s  |  time: %f ||   E_norm: %f  , dE-orth: %f ||   F=%f\n', ...
                            h, t_cur,round(scal(E1, E1), 3), round(scal(dE, E1), 3),track(end));
                end
                if sum(sqrt(w)+eps*e < 0) > 0
                    h = h/2;
                    display("!");
                else
                    break
                end
            end
            L1_E = HodgeLW_fr(B1, B2, w, e, eps);
            L0=getL0(B1, w, e, eps);
            newval=getFk_l2_con_new(L1_E, k, L0, alpha, mu);
            if  newval < track(end)
                h = h*beta; 
                break
            else
                h = h/beta;  
            end
            if h<1e-10
                h=h0;
                dE = getDotE_con_new(B1, B2, w, e, eps, k, thr, alpha, mu);
                E1 = diag(e)-h*dE;
                e = diag(E1)+nu*rand(size(e,1), 1);
                e = e/norm(e, 2);
                L1_E = HodgeLW_fr(B1, B2, w, e, eps);
                L0=getL0(B1, w, e, eps);
                newval=getFk_l2_con_new(L1_E, k, L0, alpha, mu);
                display("jump, FUCK YOU!"); display(ts(end));
                break
            end
        end
        log=[log, sort(eig(L1_E, 'vector'))];
        track=[track, newval];
        %dE = getDotE_con_new(B1, B2, w, e, eps, k, thr, alpha, mu);
        t_cur =t_cur + h;
        ts=[ts t_cur];
    end
    
    if draw
        simpleDrawB1(B1, w, e, eps)
    end

end