function [Em, Et, E_total] = bounding(K, M, h, g, D, T, gamma, N0, miu, node)
beta=0.5; % performance loss due to modulation
eta=1/1.1; % backscatter efficiency
a=2; % velocity
alpha=[0.29, 7.4]; % motion parameters

A=zeros(K,M);
B=zeros(K,1);
for k=1:K
    for m=1:M
        A(k,m)=beta*eta*norm(h(k,m))^2*norm(g(k,m))^2/N0;
    end
end

z=zeros(M,1);
z(1:length(node))=node;
C=sum(z);

if length(node)==M % leaf
    if C>=2
        % solve the TSP using Mosek
        cvx_begin quiet
        cvx_solver Mosek
        % original variables
        variable x(M,M) binary
        variable u(M)
        
        minimize trace(transpose(D)*x)
        subject to
        
        % divergence constraint
        for m=1:M
            sum(x(m,:))==z(m);
            sum(x(:,m))==z(m);
            x(m,m)==0;
        end

        % subtour elimination constraint
        for m=2:M
            for j=2:M
                if m~=j && z(m)==1 && z(j)==1
                    u(m)-u(j)+(C-1)*x(m,j)+(C-3)*x(j,m)<=C-2;
                end
            end
        end
        
        for m=2:M
            if z(m)==1
                u(m)>=1;
                u(m)<=C-1;
            else
                u(m)==0;
            end
        end
        
        cvx_end
        Tm=trace(transpose(D)*x)./a;
        Em=(alpha(1)+alpha(2)*a)*Tm;
    else
        Tm=0;
        Em=0;
    end
    
    Upsilon=T-Tm;
    z=ones(M,1);
    z(1:length(node))=node;
    
    temp=zeros(K,M);
    for k=1:K
        for m=1:M
            temp(k,m)=z(m)*A(k,m);
        end
        B(k)=max(temp(k,:));
    end
    
    [q, s] = bisection(K, gamma, B, Upsilon, miu);  
    Et=sum(s.*q);
    E_total=miu*Em+(2-miu)*Et;
else % parents
    if C>=2
        % solve the bipartite matching using Mosek
        cvx_begin quiet
        cvx_solver Mosek
        % original variables
        variable x(M,M) binary
        
        minimize trace(transpose(D)*x)
        subject to
        
        % divergence constraint
        for m=1:length(node)
            sum(x(m,:))==z(m);
            sum(x(:,m))==z(m);
            x(m,m)==0;
        end

        cvx_end
        
        Tm=trace(transpose(D)*x)./a;
        Em=(alpha(1)+alpha(2)*a)*Tm;
    else
        Tm=0;
        Em=0;
    end
    Upsilon=T-Tm;
    
    z=ones(M,1);
    z(1:length(node))=node;
    
    temp=zeros(K,M);
    for k=1:K
        for m=1:M
            temp(k,m)=z(m)*A(k,m);
        end
        B(k)=max(temp(k,:));
    end
    
    [q, s] = bisection(K, gamma, B, Upsilon, miu);  
    Et=sum(s.*q);
    E_total=miu*Em+(2-miu)*Et;
  
end

end

