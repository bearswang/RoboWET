% Load experimental data of powercast P2110 at distances of 15,14,...,3,2,1 m
% Reference:
% [1] Powercast wireless power calculator. http://www.powercastco.com/power-calculator/.
% [2] S. Wang, M. Xia, K. Huang, and Y.-C. Wu,
% “Wirelessly powered two-way communication with nonlinear energy harvesting model: Rate regions under fixed and mobile relay,”
% IEEE Transactions on Wireless Communications, 2017.
Pin =[0.0036, 0.041, 0.048, 0.056, 0.067, 0.081, 0.1,   0.127, 0.166, 0.226, 0.325, 0.507, 0.902, 2.030, 8.119]; % in mW
Pdc=[0,      0,     0,     0,     0.001, 0.011, 0.026, 0.045, 0.070, 0.107, 0.159, 0.238, 0.413, 1.097, 4.186]; % in mW
Pin_dBm=10*log10(Pin); % in dBm
Pdc_dBm=10*log10(Pdc); % in dBm
% Plot the data
hold on;
plot(Pin_dBm,Pdc_dBm,'square');

P0=0.064*1e-3; % in W [1]
Pmax=4.927*1e-3; % in W; Accordingly, the maximum input power for powercast is 10 mW [1].
nu_vec=-1:0.01:1; % searching values for nu
tau_vec=-500:1:500; % searching values for tau

% Compute the MSEs
MSE=zeros(length(nu_vec),length(tau_vec));
for i=1:length(nu_vec)
    for j=1:length(tau_vec)
        nu=nu_vec(i);
        tau=tau_vec(j);
        F1=@(x) max(Pmax./exp(-tau*P0+nu)*((1+exp(-tau*P0+nu))./(1+exp(-tau*x+nu))-1),0); % See equation (4) in Section II-D of [2]
        for l=5:length(Pin) % The term 5 is due to the first four output powers being zero
            MSE(i,j)=MSE(i,j)+abs(10*log10(F1(Pin(l)*1e-3)/1e-3)-Pdc_dBm(l))^2; % MSE with the unit of powers being dBm
        end
    end
end

% Find the minimum MSE and the corresponding index
index=find(MSE==min(min(MSE)));
[i1,j1] = ind2sub([length(nu_vec),length(tau_vec)],index);
nu=nu_vec(i1);
tau=tau_vec(j1);
F=@(x) max(Pmax*((1+exp(-tau*x+nu))^(-1)-1/(1+exp(-tau*P0+nu)))/(1-1/(1+exp(-tau*P0+nu))),0);

% Plot the model
x_dBm=-11.9:0.1:15;
x=10.^(x_dBm./10);
y=zeros(length(x),1);
for i=1:length(x)
    y(i)=F(x(i)*1e-3);
end

y_dBm=10*log10(y/1e-3);
plot(x_dBm,y_dBm,'-k');



