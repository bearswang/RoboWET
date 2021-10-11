% Define the nonlinear model
P0=0.064*1e-3;
Pmax=4.927*1e-3;
nu=0.29;
tau=274;
F=@(x) max(Pmax./exp(-tau*P0+nu)*((1+exp(-tau*P0+nu))./(1+exp(-tau*x+nu))-1),0);

chi_min=0; % Lower bound
chi_max=1; % Upper bound
ITER_MAX=1000; % Maximum iterations
chi=0.5; % Initial trial point
% Bisection search
for iter=1:ITER_MAX
    t=chi; % Trial point
    func=@(x) F(x)-t*x; % Difference between the model and the bound
    flag=0;
    for s=0:0.00001:0.1
        if func(s)>0
            flag=1; % Not a valid bound
        end
    end
     
    if flag==0 % A valid bound
        chi_max=t;
        chi=(chi_max+chi_min)./2;
    else % Not a valid bound, 
        chi_min=t; 
        chi=(chi_max+chi_min)./2;
    end
    
    if chi_max-chi<=0.001
        break;
    end
end
eta=chi;
% Defind the bound
G=@(x) min(Pmax,eta*x);

% Plot the result
x=0:0.1:30; % in mW
y1=zeros(length(x),1);
y2=zeros(length(x),1);
for i=1:length(x)
    y1(i)=F(x(i)*1e-3)*1e3; % in mW
    y2(i)=G(x(i)*1e-3)*1e3; % in mW
end
hold on;
plot(x,y2,'-r');
plot(x,y1,'-k');





