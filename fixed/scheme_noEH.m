function R = scheme_noEH(h, g, T, Er, E, Ec, Pc, Pr, Pu, miu)
% upper bound
R_ub=0.5*min(log2(1+2*(Er-Ec)/T*norm(g(:,:,2).^2)./(Pu))./miu(1),log2(1+2*(Er-Ec)/T*norm(g(:,:,1).^2)./(Pu))./miu(2)); % eqn (11)
% CF region: W. Nam [16]
kappa_min=0;
kappa_max=R_ub;
ITER_MAX=100;
kappa=zeros(ITER_MAX);
kappa(1)=kappa_max./2;

for x=1:ITER_MAX
    R=kappa(x);
    R_12=2*miu(1)*R;
    R_21=2*miu(2)*R;
    R_12_exp = 2.^(R_12);
    R_21_exp = 2.^(R_21);
    
    R_ir_exp=[R_12_exp R_21_exp]';
    R_ri_exp=[R_21_exp R_12_exp]';
    
    % optimal
    P= solve_F1_noEH(h, g, T, E, Pc, Pr, Pu, R_ir_exp, R_ri_exp);
    
    if P<=2*(Er-Ec)/T
        kappa(x+1)=(R+kappa_max)./2;
        kappa_min=R;
    else
        kappa(x+1)=(R+kappa_min)./2;
        kappa_max=R;
    end
    
    if kappa_max-kappa_min<=0.01
        break;
    end
    
    
end
R=kappa_min;
end

