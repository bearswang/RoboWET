clear all;
% Information for the platform
disp('***** This platform is used to simulate Algorithm 1 in the following paper.');
disp('Authors: S. Wang, M. Xia, K. Huang, and Y.-C. Wu');
disp('Title: Wirelessly Powered Two-Way Communication with Nonlinear Energy Harvesting Model: Rate Regions under Fixed and Mobile Relay');
disp('Publication: IEEE Trans. Wireless Commun.');

load('H.mat'); % load channels
% noises
sigmar=-70; % receiver noise in dBm
Pr=0.001*10.^(sigmar/10);
sigmau=-70; % receiver noise in dBm
Pu=0.001*10.^(sigmau/10);
sigmaz=-104; % thermal noise in dBm
Pz=0.001*10.^(sigmaz/10);
% transmission time
T=50; % in s
% local energy 
Er=300; % in J
Ec=100; % in J
E1_vec=0.4+0.4.*rand(monte,1); % in J
E2_vec=0.4+0.4.*rand(monte,1); % in J

% ciruit power
Pc = 10;% in dBm
Pc=0.001*10.^(Pc/10);
% miu vector to obtain the region
miu_vec=[0 0.2 0.4 0.5 0.6 0.8 1];
R_1=zeros(monte, length(miu_vec));
R_2=zeros(monte, length(miu_vec));

for s=1:1:length(miu_vec)
    miu(1)=miu_vec(s);
    miu(2)=1-miu_vec(s);
    tic;
    fprintf('Starting simulating miu=%d......\n',miu_vec(s));
    for m = 1: monte
        if mod(m,50)==0
            disp('50 monte carlo finished!');
            toc;
        end
        Pt=2*Er/T; % 8W
        % local power
        E=zeros(2,1);
        E(1)=E1_vec(m);
        E(2)=E2_vec(m);
        % channel
        h=h_m(:,:,:,m);
        g=g_m(:,:,:,m);
        d=d_m(:,m);

        R_1(m,s)= proposed(h, g, T, Er, E, Ec, Pc, Pr, Pu, Pz, miu);
        R_2(m,s)= scheme_noEH(h, g, T, Er, E, Ec, Pc, Pr, Pu, miu);
        
    end
end
Avg_R1 = sum(R_1,1)./monte;
Avg_R2 = sum(R_2,1)./monte;

save('region.mat');






