function R_opt = scheme_outer(N, M, h, g, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T)
% motion parameters
EM=90.54;
s=[1.5 3 1.5 0];
tM=sum(s);
% parametes in the energy harvesting model
eta=0.5732; % power conversion efficiency
Pmax=4.927*1e-3;

% channel matrix
A=zeros(N,N,2,M);
for i=1:2
    for m=1:M
        A(:,:,i,m)=g(:,:,i,m)*g(:,:,i,m)';
    end
end

% solve the outer bound problem P3 on page 22
cvx_begin sdp quiet
cvx_solver sedumi % USE Mosek for better performance!
% original variables
variable R
variable t(M)
variable V(N,N,M) hermitian
variable Q(2,M)
% slack variables
variable r(2,M)
variable X(2,M)
variable Y(2,M)

maximize R
subject to

% data rate constraint
for i=1:2
    T*miu(i)*R-sum(r(i,:))<=0; % (41a)
end

% SNR constraint
for i=1:2
    for m=1:M       
        r(i,m)+1./log(2)*rel_entr(t(m),t(m)+X(i,m))<=0; % (41b)
        r(3-i,m)+1./log(2)*rel_entr(t(m),t(m)+Y(i,m))<=0; % (41c)
    end
end

for i=1:2
    for m=1:M
    X(i,m)<=norm(h(:,:,i,m)).^2./Pr*Q(i,m);
    Y(i,m)<=real(trace(A(:,:,i,m)*V(:,:,m)))./(Pu+Pz);
    end
end

% energy constraint
for i=1:2
    for m=1:M
        temp=0;
        for l=1:m
            temp=temp+Q(i,l)+max(-eta*real(trace(A(:,:,i,l)*V(:,:,l))),-t(l)*Pmax)+(s(l)+2*t(l))*Pc;
        end
        temp<=E(i); % (41d)
    end
end
EM+real(trace(sum(V,3)))<=Er-Ec; 

tM+2*sum(t)==T; % time constraint

% basic constraints
R>=0;
for m=1:M
    t(m)>=0;
    V(:,:,m)>=0;
    Q(:,m)>=0;
    r(:,m)>=0;
end

cvx_end

if isnan(R)
    R=0;
end

R_opt = R;

end

