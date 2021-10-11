function [E_total, v0, E_iter] = algorithm2(K, M, h, g, D, T, gamma, N0, miu, L, ITER_MAX)
beta=0.5; % performance loss due to modulation
eta=1/1.1; % backscatter efficiency
a=2; % velocity in m/s
alpha=[0.29, 7.4]; % UGV parameters

A=zeros(K,M); 
for k=1:K
    for m=1:M
        A(k,m)=beta*eta*norm(h(k,m))^2*norm(g(k,m))^2/N0; % equation (7)
    end
end

B=zeros(K,1); % B in equation (11)

E_iter=zeros(ITER_MAX,1); % result matrix
v0=zeros(M,1); 
v0(1)=1; % no UGV movement
v=v0; % initial v

% start iterative procedure of randomized local search
for iter =1:ITER_MAX
    C=sum(v); % number of places to visit
    if C==1 % no UGV movement
        Tm=0; % 0 moving time
        Em=0; % 0 moving energy
    else
        J=1e6;
        % solve problem (13) using Mosek
        cvx_begin quiet
        cvx_solver Mosek
        variable x(M,M) binary
        variable u(M)
        
        minimize trace(transpose(D)*x) % minimize path length
        subject to
        
        % divergence constraint
        for m=1:M
            sum(x(m,:))==v(m); % (6c)
            sum(x(:,m))==v(m); % (6c)
            x(m,m)==0;
        end
        
        % subtour elimination constraint
        for m=2:M
            for j=2:M
                if m~=j
                    u(m)-u(j)+(C-1)*x(m,j)+(C-3)*x(j,m)<=C-2+J*(2-v(m)-v(j)); % (6d)
                end
            end
        end
        
        for m=2:M
            if v(m)==1
                u(m)>=1; % (6e)
                u(m)<=C-1; % (6e)
            else
                u(m)==0;
            end
        end
        
        cvx_end
        Tm=trace(transpose(D)*x)./a; % moving time
        Em=(alpha(1)+alpha(2)*a)*Tm; % moving energy
    end
    
    Upsilon=T-Tm; % equation (14)

    temp=zeros(K,M);
    for k=1:K
        for m=1:M
            temp(k,m)=v(m)*A(k,m);
        end
        B(k)=max(temp(k,:));
    end
    [q, s] = bisection(K, gamma, B, Upsilon, miu);

    Et=sum(q.*s);
    E_current=miu*Em+(2-miu)*Et;
    
    if iter==1
        E_iter(iter)=E_current;
    else
        if E_current<=E_iter(iter-1)
            v0=v;
            E_iter(iter)=E_current;
        else
            E_iter(iter)=E_iter(iter-1);
        end

    end
    
    while 1
        v=v0;
        for coin=1:L
            temp=1+ceil((M-1)*rand);
            v(temp)=1-v0(temp);
        end
        if sum(v)<=K+1
            break;
        end
    end

    E_total=min(E_iter);
    
end

