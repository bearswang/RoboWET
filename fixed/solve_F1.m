function P_opt = solve_F1(h, g, T, E, Pc, Pr, Pu, Pz, R_ir_exp, R_ri_exp)
u=h(:,:,2)-h(:,:,1)'*h(:,:,2)./norm(h(:,:,1)).^2*h(:,:,1); % u is defined under equation (18)
u=u./norm(u); % normalize u

temp=zeros(2,1);
temp(1)=R_ir_exp(1)./(R_ir_exp(1)+R_ir_exp(2));
temp(2)=R_ir_exp(2)./(R_ir_exp(1)+R_ir_exp(2));
alpha=R_ir_exp-temp; % alpha is defined above equation (13)

gamma_vec=0:0.05:1;
p=zeros(length(gamma_vec),1);
for s=1:length(gamma_vec)
    gamma=gamma_vec(s);
    p(s) = fgamma(h, g, T, E, Pc, Pr, Pu, Pz, R_ri_exp, gamma, alpha, u);
end
P_opt=min(p);

end

