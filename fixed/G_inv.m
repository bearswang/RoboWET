function y = G_inv(x)
P0=0.064*1e-3;
Pmax=4.927*1e-3;
nu=0.29;
tau=274;

if x<=0
    y=0;
else
    if x>=Pmax
        y=10^10;
    else
        temp=exp(-tau*P0+nu);
        y=tau^(-1)*(nu-log((1+temp)./(1+Pmax^(-1)*temp*x)-1));
    end
    
end

