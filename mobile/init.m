function [R0,t0,w0,a0,c0,Q0,feasible]= init(N, M, h, g, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T)
% init for DC programming
EM=90.54;
s=[1.5 3 1.5 0];
tM=sum(s);
eta=0.5732; % power conversion efficiency
Pmax=4.927*1e-3;
tau=274;
A=zeros(N,N,2,M);
for i=1:2
    for m=1:M
        A(:,:,i,m)=g(:,:,i,m)*g(:,:,i,m)';
    end
end

% solve the outer bound problem P3 on page 22
cvx_begin sdp quiet
cvx_solver sedumi
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
    feasible = 0;
    V0=zeros(N,N,M);
else
    feasible = 1;
    V0=V;
end


t0=(T-sum(s))/2/M*ones(M,1); % equal initial time allocation

w0=zeros(N,1,M); % intial receiver
for m=1:M
    z=h(:,:,2,m)-h(:,:,1,m)'*h(:,:,2,m)./norm(h(:,:,1,m)).^2*h(:,:,1,m); % t is defined under eqn (19)
    z=z./norm(z); % normalize t
    w0(:,:,m)=sqrt(0.5)*h(:,:,1,m)./norm(h(:,:,1,m))+sqrt(0.5)*exp(sqrt(-1)*angle(h(:,:,2,m)'*h(:,:,1,m)))*z; % eqn (38)
    w0(:,:,m)=w0(:,:,m)./norm(w0(:,:,m)); % the receiver
end

omega0=zeros(2,M);
for i=1:2
    for m=1:M
        omega0(i,m)=0.1*real(trace(A(:,:,i,m)*V0(:,:,m)));
    end
end

a0=zeros(2,M);
for i=1:2
    for m=1:M
        a0(i,m)=(Pu./real(trace(A(:,:,i,m)*V0(:,:,m)))+Pz./omega0(i,m))^(-1);
    end
end

c0=zeros(2,M);
for i=1:2
    for m=1:M
        c0(i,m)=t0(m)*exp((real(trace(A(:,:,i,m)*V0(:,:,m)))-omega0(i,m))*tau/t0(m));
    end
end

snr=100; % 20 dB initial uplink SNR
Q0=zeros(2,M);
for m=1:M
    Q0(:,m)=snr*t0(m)*Pr;
end
R0=0;
end

