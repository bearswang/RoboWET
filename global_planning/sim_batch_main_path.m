clear all;
% Information for the platform
disp('***** This platform is used to simulate the algorithm in the following paper.');
disp('S. Wang, M. Xia, Y.-C. Wu');
disp('Backscatter data collection with unmanned ground vehicle');
disp('IEEE Transactions on Wireless Communications');
load('map.mat')

% time budget
T=500; % in s
% data rate requirement
miu=1;
sigma=-90;
N0=0.001*10.^(sigma/10);
gamma=1+rand(K,1);

v_all=zeros(M,3,monte);
mon=1; % index of map realization
fprintf('Starting simulating path......');
% load channels
h=h_all(:,:,mon);
g=g_all(:,:,mon);
D=D_all(:,:,mon);

[E_opt, Em_opt, Et_opt, v_opt, history, E_ls, v_ls, E_iter] = proposed(K, M, h, g, D, T, gamma, N0, miu, 20);
v_all(:,1,mon)=ones(M,1); % full path
v_all(:,2,mon)=v_ls; % path of Algorithm 2
v_all(:,3,mon)=v_opt; % optimal path obtained from Algorithm 1
% store the results
txt = sprintf('example_%d.mat',sigma);
save(txt);


