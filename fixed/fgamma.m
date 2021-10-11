function P_opt = fgamma(h, g, T, E, Pc, Pr, Pu, Pz, R_ri_exp, gamma, alpha, u)
theta=R_ri_exp-1;
w=sqrt(gamma)*h(:,:,1)./norm(h(:,:,1))+sqrt(1-gamma)*exp(sqrt(-1)*angle(h(:,:,2)'*h(:,:,1)))*u; % eqn (18)
if norm(w)==0
   P_opt=10^10;
   return;
end
w=w./norm(w); % normalize the receiver

B=zeros(2,1);
for i=1:2
    B(i)=G_inv((alpha(i)*Pr)./norm(w'*h(:,:,i))^2+2*Pc-2*E(i)/T);
end
% derive the optimal power splitter
splitter=zeros(2,1); % power splitter
for i=1:2
if theta(i)==0
    splitter(i)=0;
else
    if B(i)==0
        splitter(i)=1;
    else
        temp=Pz-Pu+B(i)./theta(i);
        splitter(i)=2*Pz./(temp+sqrt(temp^2+4*Pu*Pz));
    end
end
end
% derive the optimal transmitter
xi=zeros(2,1); % xi defined in eqn (24)
for i=1:2
   if theta(i)==0
       xi(i)=B(i);
   else
       if B(i)==0
           xi(i)=theta(i)*(Pu+Pz);
       else
           xi(i)=theta(i)*(Pu+Pz./splitter(i));
       end
   end    
end
z=g(:,:,2)-g(:,:,1)'*g(:,:,2)./norm(g(:,:,1)).^2*g(:,:,1); % z defined above (16)
A=zeros(3,2); 
p=zeros(3,1);
rho=zeros(3,1); % there are three candidates, See eqn (27)
rho(1)=1; % See eqn (27)
rho(2)=norm(g(:,:,2)'*g(:,:,1))^2./(norm(g(:,:,1))^2*norm(g(:,:,2))^2); % See eqn (27)
denominator=(sqrt(xi(2))*norm(g(:,:,1))-sqrt(xi(1))*norm(g(:,:,2)'*g(:,:,1))./norm(g(:,:,1)))^2+xi(1)*norm(z)^2; % See eqn (28)
rho(3)=xi(1)*norm(z)^2./denominator; % intersection point, see eqn (28)

for j=1:3
    A(j,1)=rho(j)*norm(g(:,:,1))^2; % A defined in (20)
    A(j,2)=(sqrt(rho(j))*norm(g(:,:,2)'*g(:,:,1))./norm(g(:,:,1))+sqrt(1-rho(j))*norm(z))^2; % A defined in (20)
    p(j)=max(xi(1)./A(j,1),xi(2)/A(j,2));
end 

P_opt=min(p); % the optimal relay power
     
end

