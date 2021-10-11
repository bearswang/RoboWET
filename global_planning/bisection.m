function [q, s] = bisection(K, gamma, B, Upsilon, miu)
if Upsilon<=0
    s=1;
    return;
end

s=zeros(K,1); % transmission time allocation
Lambda_inverse=@(x) (1+log(2)*2.^(1./x)./x-2.^(1./x)); % -\nabla\Theta, i.e., inverse function of Lambda in (16)
term=Lambda_inverse(Upsilon./sum(gamma));
rho_min=(2-miu)./sum(gamma)./(sum(gamma.*B))*term;
rho_max=(2-miu)./min(B)*term;

ITER_MAX=1000; % Maximum iterations
rho=(rho_min+rho_max)/2; % Initial trial point

for iter=1:ITER_MAX
    
    % bisection to solve equation (15) given fixed rho
    s_min=zeros(K,1);
    s_max=Upsilon*ones(K,1);
    s_try=(s_min+s_max)/2;
    
    for ii=1:ITER_MAX
        rho_hat=(2-miu)*Lambda_inverse(s_try./gamma);
        diff=rho_hat-rho;
        s_min= (diff>=0).*(s_try-s_min)+s_min;
        s_max= (diff<0).*(s_try-s_max)+s_max;
        s_try=(s_max+s_min)./2;
        if norm(s_max-s_try)<=1e-6
            break;
        end
    end
    s=s_try;

    if sum(s)>=Upsilon
        rho_min=rho;
        rho=(rho_max+rho_min)./2;
    else
        rho_max=rho;
        rho=(rho_max+rho_min)./2;
    end
    if norm(rho_max-rho)<=1e-6
        break;
    end
    
end

q=1./B.*(2.^(gamma./s)-1);

end