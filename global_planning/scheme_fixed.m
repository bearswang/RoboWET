function [Em,Et,Etotal]= scheme_fixed(K, M, h0, g0, T, gamma, N0, miu)  
beta=0.5;
eta=1/1.1;
A=zeros(K,1);
for k=1:K
    A(k)=beta*eta*norm(h0(k))^2*norm(g0(k))^2/N0;
end

[p, t] = bisection(K, gamma, A, T, miu);
Et=sum(t.*p);

Em=0;
Etotal=miu*Em+(2-miu)*Et;

end

