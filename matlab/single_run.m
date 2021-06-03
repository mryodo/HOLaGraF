function [ans, track, ts, ending, history, history_de, log, st_p, G]=single_run(G, k, thrs, h0, my_beta, thr0, fileID)

    h=h0; fl=0;
    L1 = HodgeLW_fr(G);
    p=sum(abs(L1) * ones(size(L1, 1), 1)<1e-3);

    track = [getFk_l2(G, k, p, thrs, fileID)];
%     if getFk2(G, k, p, thrs, fileID)>1e-6
%         cs = linspace(0, 0.25, 10);
%         dE=getGradNum(G, k, p, thrs, thr0, 1, fileID);
%         cor=0;
%         for i=1:10
%            G1=G; G1.e=G.e-cs(i)*dE/norm(dE, 2);
%            G1.e=G1.e/norm(G1.e, 2);
%            if getFk2(G, k, p, thrs, fileID)<1e-6
%               fprintf("Correction with: %f", cs(i));
%               G.e=G1.e;
%               cor=1;
%               break
%            end
%         end
%         if ~cor
%             %fprintf("alph: %f \n Penalization from the start! Abort, abort!\n", thrs.alph);
%             %return         
%         end
%     end
    t_cur = 0; ts = [0];
    log = []; h_log=[]; h_desc=[];
    history=[G.e]; history_de=[];

    ending=0; step_num=0; st_p=0; jump=0;

    while (track(end)>1e-5) && (fl==0) && (step_num<=1000)
         e0 = G.e;
         dE=getGrad(G, k, p, thrs, thr0, 1, fileID);
         history_de=[history_de dE];
         if (max(abs(dE))<1e-4)
            fl=1;
            ending=1;
            break
         end
         
         while 1
             e1 = G.e-h*dE;
             e1(sqrt(G.w)+G.eps0*e1<0)=-1./G.eps0*sqrt(G.w(sqrt(G.w)+G.eps0*e1<0));
             e1=e1/norm(e1, 2);
             G1=G; G1.e=e1;
             newval=getFk_l2(G1, k, p, thrs, fileID);
             L1_E=HodgeLW_fr(G1);
%              if p ~= sum(abs(L1_E) * ones(size(L1_E, 1), 1)<1e-3)
%                  p=sum(abs(L1_E) * ones(size(L1_E, 1), 1)<1e-3);
%                  st_p=st_p+1;
%                  newval=getFk_l2(G1, k, p, thrs, fileID);
%                  G1.e=e1;
%                  jump=1;
%                  break
%              end
             if  (newval < track(end)) || (jump==1)      
                 if jump==0
                     h=h*my_beta;
                 end
                 h_log=[h_log h]; h_desc=[h_desc 1];
                 G.e = G1.e;
                 jump=0;
                 break
             else
                 G.e=e0;
                 h=h/my_beta;
                 h_log=[h_log h]; h_desc=[h_desc 2];
             end
             if (h<1e-10)
                fl=1; ending=2;
                break
             end
         end
         
         log=[log sort(eig(L1_E))];
         track=[track newval];
         history=[history G.e];
         t_cur = t_cur + h; ts=[ts t_cur];
         step_num = step_num+1;
     end
     ans=history(:, end);
end