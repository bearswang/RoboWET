function [E_opt, Em_opt, Et_opt, v_opt, history, E_ls, v_ls, E_iter, p_max] = proposed(K, M, h, g, D, T, gamma, N0, miu, ITER_ls)
L=3; % step-size for Algorithm 2
[E_ls, v_ls, E_iter]= algorithm2(K, M, h, g, D, T, gamma, N0, miu, L, ITER_ls); % initialization by running Algorithm 2
[E_opt, v_opt, history, Em_opt, Et_opt]= algorithm1(K, M, h, g, D, T, gamma, N0, miu, E_ls, ITER_ls); % branch and bound
p_max = peak_power(K, M, h, g, D, T, gamma, N0, miu, v_opt); % peak power
end

