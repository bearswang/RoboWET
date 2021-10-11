clear all;
% Information for the platform
disp('***** This platform is used to simulate the algorithm in the following paper.');
disp('S. Wang, M. Xia, Y.-C. Wu');
disp('Backscatter data collection via unmanned ground vehicle');
disp('IEEE Transactions on Wireless Communications');
load('map.mat'); % load the map and channels

T=500; % time budget, in second
miu=1; % weighting factor to balance motion versus communication
sigma_vec=[-140, -120, -100, -80, -60]; % noise power in dBm
gamma=1+rand(K,1); % data collection target in bit/Hz
ITER_MAX=50; % number of iterations for local search

% energy for the fixed scheme
Em1=zeros(monte, length(sigma_vec));
Et1=zeros(monte, length(sigma_vec));
E1=zeros(monte, length(sigma_vec));
% energy for the full path scheme
Em2=zeros(monte, length(sigma_vec));
Et2=zeros(monte, length(sigma_vec));
E2=zeros(monte, length(sigma_vec));
% energy for the proposed algorithm
Em3=zeros(monte, length(sigma_vec));
Et3=zeros(monte, length(sigma_vec));
E3=zeros(monte, length(sigma_vec));

Pmax=zeros(monte, length(sigma_vec)); % peak power

for s=1:1:length(sigma_vec)
    sigma=sigma_vec(s); % noise
    N0=0.001*10.^(sigma/10); % noise power
    tic;
    fprintf('Starting simulating noise=%d......\n',sigma_vec(s));
    for mon = 1: monte
        if mod(mon,50)==0
            disp('50 monte carlo finished!');
            toc;
        end
        
        % load channel at mon
        h=h_all(:,:,mon);
        g=g_all(:,:,mon);
        D=D_all(:,:,mon);
        % channels at the starting point
        h0=h(:,1);
        g0=g(:,1);     
        
        [Em1(mon,s), Et1(mon,s), E1(mon,s)]= scheme_fixed(K, M, h0, g0, T, gamma, N0, miu); % no UGV movement
        [Em2(mon,s), Et2(mon,s), E2(mon,s)]= scheme_fullpath(K, M, h, g, D, T, gamma, N0, miu); % full path visiting all vertices
        % proposed Algorithm 1
        [E3(mon,s), Em3(mon,s), Et3(mon,s), v_opt, history, E_ls, v_ls, E_iter, p_max] = proposed(K, M, h, g, D, T, gamma, N0, miu, 20); 

        Pmax(mon,s)=p_max; % track the largest transmit power
        
    end
end

Avg_E1 = sum(E1,1)./monte;
Avg_E2 = sum(E2,1)./monte;
Avg_E3 = sum(E3,1)./monte;
Pmax_all = max(Pmax);

save('noise.mat');






