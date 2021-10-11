function R_iter= scheme_inner(N, M, h, g, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T, ITER_MAX)
% moving parameters
EM=90.54;
s=[1.5 3 1.5 0];
% parameters in the energy harvesting model
P0=0.064*1e-3;
Pmax=4.927*1e-3;
nu=0.29;
tau=274;
% terms involved in function H
H_term1=Pmax/exp(-tau*P0+nu);
H_term2=1+exp(-tau*P0+nu);

A=zeros(N,N,2,M);
for i=1:2
    for m=1:M
        A(:,:,i,m)=g(:,:,i,m)*g(:,:,i,m)';
    end
end

% initialization from problem P2[0]
[R0,t0,w0,a0,c0,Q0,feasible]=feasible_point_pursuit(N, M, h, g, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T);

R_iter=zeros(ITER_MAX+1,1); % metric to track the iteration
if feasible==0 % not feasible
    return;
else
    R_iter(1)=R0;
end

% solving problem P2[n+1]
for iter=1:ITER_MAX
    % the coefficients in (34)
    Phi_term1=zeros(M,1);
    Phi_term2=zeros(M,1);
    Phi_term3=zeros(M,1);
    for m=1:M
        Phi_term1(m)=-t0(m)*log2(1+t0(m)*Pr./(Q0(1,m)+Q0(2,m)));
        Phi_term2(m)=1./log(2)*t0(m)^2./(t0(m)*(Q0(1,m)+Q0(2,m))+(Q0(1,m)+Q0(2,m))^2./Pr);
        Phi_term3(m)=-log2(1+t0(m)*Pr./(Q0(1,m)+Q0(2,m)))-t0(m)/log(2)./(t0(m)+(Q0(1,m)+Q0(2,m))./Pr);
    end
    
    Psi_term1=zeros(2,M);
    Psi_term2=zeros(2,M);
    Psi_term3=zeros(2,M);
    for i=1:2
        for m=1:M
            Psi_term1(i,m)=t0(m)/tau*log(c0(i,m)/t0(m));
            Psi_term2(i,m)=t0(m)./(tau*c0(i,m));
            Psi_term3(i,m)=1/tau*log(c0(i,m)/t0(m))-1/tau;
        end
    end
    
    % solve the inner bound problem P2[0] on page 19
    cvx_begin sdp quiet
    cvx_solver sedumi
    % original variables
    variable R
    variable t(M)
    variable V(N,N,M) hermitian
    variable Q(2,M)
    variable w(N,1,M) complex
    variable omega(2,M)
    
    % slack variables
    variable r(2,M)
    variable a(2,M)
    variable b(2,M)
    variable c(2,M)
    variable H_term3(2,M)
    
    maximize R
    subject to
    
    for i=1:2
        T*miu(i)*R-sum(r(i,:))<=0; % (line 1)
    end
    
    for i=1:2
        for m=1:M
            Phi_tilde=Phi_term1(m)+Phi_term2(m)*(Q(1,m)+Q(2,m)-Q0(1,m)-Q0(2,m))+Phi_term3(m)*(t(m)-t0(m));
            r(i,m)+1./log(2)*rel_entr(t(m),Q(i,m)/Pr)+Phi_tilde<=0; % (line 2)
        end
    end
    
    for i=1:2
        for m=1:M
            r(3-i,m)+1./log(2)*rel_entr(t(m),t(m)+a(i,m))<=0; % (line 3)
        end
    end
    
    for i=1:2
        for m=1:M
            Upsilon_tilde=-2./a0(i,m)+1./a0(i,m).^2*a(i,m);
            Pu*pow_p(real(trace(A(:,:,i,m)*V(:,:,m))),-1)+Pz*pow_p(omega(i,m),-1)+Upsilon_tilde<=0; % (line 4)
        end
    end
    
    for i=1:2
        for m=1:M
            temp=0;
            for l=1:m
                H=H_term1*t(l)-H_term1*H_term2*t(l)+H_term1*H_term2*exp(nu)*H_term3(i,l);
                temp=temp+b(i,l)+H+(2*t(l)+s(l))*Pc;
            end
            temp<=E(i); % (line 5)
        end
    end
    
    for i=1:2
        for m=1:M
            [H_term3(i,m) t(m);t(m) c(i,m)+exp(nu)*t(m)]>=0; % (line 6)
        end
    end
    
    for i=1:2
        for m=1:M
                Xi_tilde=-2*real(w0(:,:,m)'*h(:,:,i,m)*h(:,:,i,m)'*w(:,:,m))/Q0(i,m)+real(w0(:,:,m)'*h(:,:,i,m)*h(:,:,i,m)'*w0(:,:,m))/(Q0(i,m).^2)*Q(i,m);
                pow_p(b(i,m),-1)+Xi_tilde<=0; % (line 7)
        end
    end
    
    
    for i=1:2
        for m=1:M
            Psi_tilde=Psi_term1(i,m)+Psi_term2(i,m)*(c(i,m)-c0(i,m))+Psi_term3(i,m)*(t(m)-t0(m));
            real(trace(A(:,:,i,m)*V(:,:,m)))-omega(i,m)>=Psi_tilde; % (line 8)
        end
    end
    
    % (line 9 (g)(h))
    EM+real(trace(sum(V,3)))<=Er-Ec;
    2*sum(t)+sum(s)==T;
    % basic constraints
    R>=0;
    for m=1:M
        t(m)>=1e-3;
        V(:,:,m)>=0;
        norm(w(:,:,m))<=1;
        Q(:,m)>=0;
        c(:,m)>=t(m); 
        a(:,m)>=0;
        b(:,m)>=0;
        r(:,m)>=0;
    end
    
    cvx_end
    if ~isnan(R) % check CVX status
        t0=t;
        w0=w;
        Q0=Q;
        a0=a;
        c0=c;
        R_iter(iter+1)=R;
    else % CVX numerical error
        R_iter(iter+1:end)= R_iter(iter);
        break;
    end
    
    if norm(R_iter(iter+1)-R_iter(iter))<=0.1 % check convergence
       R_iter(iter+1:end)= R_iter(iter+1);
       break; 
    end
    
end


end

