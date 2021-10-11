function P_opt = fgamma_noEH(h, g, T, E, Pc, Pr, Pu, R_ri_exp, gamma, alpha, u)
theta=R_ri_exp-1;
w=sqrt(gamma)*h(:,:,1)./norm(h(:,:,1))+sqrt(1-gamma)*exp(sqrt(-1)*angle(h(:,:,2)'*h(:,:,1)))*u;
if norm(w)==0
   P_opt=10^10;
   return;
end
w=w./norm(w);

B=zeros(2,1);
for i=1:2
    B(i)=max(0,(alpha(i)*Pr)./norm(w'*h(:,:,i))^2+2*Pc-2*E(i)/T);
end

if ~isempty(nonzeros(B))
    P_opt=1e10;
    return;
end

xi=zeros(2,1);
for i=1:2
    xi(i)=theta(i)*Pu;
end

z=g(:,:,2)-g(:,:,1)'*g(:,:,2)./norm(g(:,:,1)).^2*g(:,:,1);
A=zeros(3,2); 
p=zeros(3,1);
rho=zeros(3,1);
rho(1)=1;
rho(2)=norm(g(:,:,2)'*g(:,:,1))^2./(norm(g(:,:,1))^2*norm(g(:,:,2))^2);
denominator=(xi(2)^0.5*norm(g(:,:,1))-xi(1)^0.5*norm(g(:,:,2)'*g(:,:,1))./norm(g(:,:,1)))^2+xi(1)*norm(z)^2;
rho(3)=xi(1)*norm(z)^2./denominator;

for i=1:3
    A(i,1)=rho(i)*norm(g(:,:,1))^2;
    A(i,2)=(sqrt(rho(i))*norm(g(:,:,2)'*g(:,:,1))./norm(g(:,:,1))+sqrt(1-rho(i))*norm(z))^2;
    p(i)=max(xi(1)./A(i,1),xi(2)/A(i,2));
end 

P_opt=min(p);
    
end

