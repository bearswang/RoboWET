function [Em,Et,Etotal] = scheme_fullpath(K, M, h, g, D, T, gamma, N0, miu)
beta=0.5;
eta=1/1.1;
alpha=[0.29, 7.4];
a=2; % 2m/s

A=zeros(K,M);
B=zeros(K,1);
for k=1:K
    for m=1:M
        A(k,m)=beta*eta*norm(h(k,m))^2*norm(g(k,m))^2/N0;
    end
    B(k)=max(A(k,:));
end
v=ones(M,1);
% solve the TSP using Gurobi
cvx_begin quiet
cvx_solver Mosek
% original variables
variable x(M,M) binary
variable u(M)

minimize trace(transpose(D)*x) % total path length
subject to

% divergence constraint
for m=1:M
    if v(m)==1;
        sum(x(m,:))==1;
        sum(x(:,m))==1;
    else
        sum(x(m,:))==0;
        sum(x(:,m))==0;
    end
    x(m,m)==0;
end
% subtour elimination constraint
for m=2:M
    for j=2:M
        if m~=j && v(m)==1 && v(j)==1
            u(m)-u(j)+(M-1)*x(m,j)+(M-3)*x(j,m)<=M-2;
        end
    end
end

for m=2:M
    if v(m)==1
        u(m)>=1;
        u(m)<=M-1;
    else
        u(m)==0;
    end
end

cvx_end


Sm=trace(transpose(D)*x)./a;
Em=(alpha(1)+alpha(2)*a)*Sm;
[p, t] = bisection(K, gamma, B, T-Sm, miu);
Et=sum(t.*p);
Etotal=miu*Em+(2-miu)*Et;

end

