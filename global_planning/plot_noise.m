load('noise.mat');
semilogy(sigma_vec,Avg_E1,'-b');
hold on;
semilogy(sigma_vec,Avg_E2,'-r');
semilogy(sigma_vec,Avg_E3,'-square');