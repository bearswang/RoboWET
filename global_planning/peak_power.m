function p_max = peak_power(K, M, h, g, D, T, gamma, N0, miu, v)
beta=0.5;
eta=1/1.1;
a=2;
alpha=[0.29, 7.4];
A=zeros(K,M);
temp=zeros(K,M);
B=zeros(K,1);
for k=1:K
    for m=1:M
        A(k,m)=beta*eta*norm(h(k,m))^2*norm(g(k,m))^2/N0;
    end
end

C=sum(v);

if C>=2
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
        sum(x(m,:))==v(m);
        sum(x(:,m))==v(m);
        x(m,m)==0;
    end
    
    % subtour elimination constraint
    % subtour elimination constraint
    for m=2:M
        for j=2:M
            if m~=j && v(m)==1 && v(j)==1
                u(m)-u(j)+(C-1)*x(m,j)+(C-3)*x(j,m)<=C-2;
            end
        end
    end
    
    for m=2:M
        if v(m)==1
            u(m)>=1;
            u(m)<=C-1;
        else
            u(m)==0;
        end
    end
    
    cvx_end
    
    Sm=trace(transpose(D)*x)./a;
else
    Sm=0;
end

for k=1:K
    for m=1:M
        temp(k,m)=v(m)*A(k,m);
    end
    B(k)=max(temp(k,:));
end

[p, ~] = bisection(K, gamma, B, T-Sm, miu);
p_max=max(p);



end

