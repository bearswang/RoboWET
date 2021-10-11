clear all;
% Information for the platform
disp('***** This platform is used to simulate the algorithm in the following paper.');
disp('S. Wang, M. Xia, Y.-C. Wu');
disp('Wirelessly powered two-way communication with nonlinear energy harvesting model: Rate regions under fixed and mobile relay');
disp('IEEE Transactions on Wireless Communications');

load('channel.mat'); % load channels
monte=1;
% noises
sigmar=-70; % in dBm
Pr=0.001*10.^(sigmar/10);
sigmau=-70; % in dBm
Pu=0.001*10.^(sigmau/10);
sigmaz=-104; %  in dBm
Pz=0.001*10.^(sigmaz/10);
% transmission time
T=50; % in s
% local energy
Er=300; % in J
Ec=100;
E1=0.8.*ones(monte,1); % in mJ
E2=0.8.*ones(monte,1); % in mJ

% ciruit power
Pc = 10;% in dBm
Pc=0.001*10.^(Pc/10);
% miu vector to obtain the region
miu_vec=[0 0.3 0.5 0.7 1];

% registers
R1=zeros(monte, length(miu_vec));
R2=zeros(monte, length(miu_vec));
R3=zeros(monte, length(miu_vec));
R4=zeros(monte, length(miu_vec));
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
        % local energy
        E=[E1(m) E2(m)]';
        
        ITER_MAX=20;
        % trajectory 1
        h1=h_m_1(:,:,:,:,m);
        g1=g_m_1(:,:,:,:,m);
        
        R = scheme_inner(N, M, h1, g1, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T, ITER_MAX);
        R1(m,s)=max(R);
        R2(m,s) = scheme_outer(N, M, h1, g1, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T);
        
        
        % trajectory 2
        h2=h_m_2(:,:,:,:,m);
        g2=g_m_2(:,:,:,:,m);
        
        R = scheme_inner(N, M, h2, g2, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T, ITER_MAX);
        R3(m,s)=max(R);
        R4(m,s) = scheme_outer(N, M, h2, g2, Er, Ec, E, Pc, Pr, Pu, Pz, miu, T);
    end
    
end

Avg_R1 = sum(R1,1)./monte;
Avg_R2 = sum(R2,1)./monte;
Avg_R3 = sum(R3,1)./monte;
Avg_R4 = sum(R4,1)./monte;

save('mobile.mat');






