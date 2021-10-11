load('mobile.mat');
hold on;
plot(Avg_R2.*miu_vec,Avg_R2.*(1-miu_vec),'-.ro');
plot(Avg_R4.*miu_vec,Avg_R4.*(1-miu_vec),'-.ko');
plot(Avg_R1.*miu_vec,Avg_R1.*(1-miu_vec),'-rsquare');
plot(Avg_R3.*miu_vec,Avg_R3.*(1-miu_vec),'-ksquare');

