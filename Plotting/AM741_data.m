%%% RMSE Data
%% Reproduced Advection Diffusion
u_n0 = [0, 1.22e-2, 8.08e-4, 4.09e-4];
u_n1 = [1.01e-4,1.22e-2,8.30e-4,5.06e-4];
u_n5 = [2.50e-3,1.41e-2,1.08e-3,5.9e-4];
u_n10 = [1.01e-2,1.65e-2,2.23e-3,1.13e-3];
u_n25 = [6.24e-2,1.88e-2,7.23e-3,7.94e-3];
u_n50 = [2.38e-1,2.67e-2,3.16e-2,5.86e-2];
u = [u_n0;u_n1;u_n5;u_n10;u_n25;u_n50];

ut_n0 = [5.39e-5,1.12,1.41e2,5.19e-1];
ut_n1 = [4.93e2,4.05,2.84e2,7.01e-1];
ut_n5 = [4.02e4,3.14e2,4.87e3,9.44e-1];
ut_n10 = [6.72e4,9.76e2,1.16e4,3.30e-1];
ut_n25 = [1.39e5,6.32e2,8.30e4,1.19];
ut_n50 = [1.97e5,9.13e3,4.24e5,9.15];
ut = [ut_n0;ut_n1;ut_n5;ut_n10;ut_n25;ut_n50];

ux_n0 = [5.77e-4,3.33e-2,2.83,3.79e-2];
ux_n1 = [3.02e-1,4.48e-2,3.24,1.6e-2];
ux_n5 = [1.12e1,4.21e-1,6.97,4.71e-2];
ux_n10 = [7.21e1,1.171,1.72e1,1e-2];
ux_n25 = [3.57e2,4.54,8.30e1,4.13e-2];
ux_n50 = [1.36e3,8.04e1,4.54e2,3.22e-1];
ux = [ux_n0;ux_n1;ux_n5;ux_n10;ux_n25;ux_n50];

uxx_n0 = [3.69e-2,4.86e1,5.56e1,4.08e-1];
uxx_n1 = [2.53e1,4.86e1,5.85e1,6.90e-1];
uxx_n5 = [9.89e2,6.05e1,1.16e2,5.39e-1];
uxx_n10 = [5.46e3,1.38e2,5.30e2,6.04e-1];
uxx_n25 = [3.75e4,1.63e2,2.59e3,6.92e-1];
uxx_n50 = [1.19e5,4.15e2,1.57e4,1.14];
uxx = [uxx_n0;uxx_n1;uxx_n5;uxx_n10;uxx_n25;uxx_n50];

plot_RMSE(u,ut,ux,uxx,'\textbf{Advection Diffusion}');

%% Reproduced Fishers

f_u_n0 = [0, 6.28e-5,3.80e-6,4.67e-4];
f_u_n1 = [9.95e-5,6.90e-5,9.95e-6,3.39e-4];
f_u_n5 = [2.49e-3,2.12e-4,1.60e-4,4.45e-4];
f_u_n10 = [1.01e-2,7.45e-4,6.57e-4,7.77e-4];
f_u_n25 = [6.23e-2,3.90e-3,4.40e-3,5.44e-3];
f_u_n50 = [2.42e-1,1.57e-2,3.15e-2,7.65e-2];
f_u = [f_u_n0;f_u_n1;f_u_n5;f_u_n10;f_u_n25;f_u_n50];

f_ut_n0 = [3.16e-5,2.83e-4,9.43e-2,6.45e-2];
f_ut_n1 = [3.97,1.63e-1,1.04e1,1.11e-2];
f_ut_n5 = [1.04e2,4.06,1.34e2,4.63e-2];
f_ut_n10 = [3.76e2,2.93e1,4.73e3,8.3e-2];
f_ut_n25 = [2.39e3,7.38e1,2.40e4,3.65e-2];
f_ut_n50 = [9.82e3,6.16e2,4.02e4,9.97e-1];
f_ut = [f_ut_n0;f_ut_n1;f_ut_n5;f_ut_n10;f_ut_n25;f_ut_n50];

f_ux_n0 = [4.46e-4,2.83e-3,5.54e-1,4.47e-2];
f_ux_n1 = [6.76e1,3.89,1.96e1,1.64e-2];
f_ux_n5 = [1.86e3,5.31e1,1.74e2,5.52e-2];
f_ux_n10 = [4.85e3,1.95e2,6.96e3,1.08e-1];
f_ux_n25 = [1.88e4,1.74e3,1.17e4,1.96e-1];
f_ux_n50 = [1.29e5,4.40e3,1.76e4,1.25];
f_ux = [f_ux_n0;f_ux_n1;f_ux_n5;f_ux_n10;f_ux_n25;f_ux_n50];

f_uxx_n0 = [4.67e-3,2.01e-1,5.17e-1,4.24e-1];
f_uxx_n1 = [3.47e3,2.52e1,3.23e1,2.51];
f_uxx_n5 = [3.37e4,3.56e1,2.25e2,4.28e-1];
f_uxx_n10 = [6.17e5,1.09e3,3.14e3,4.18e-2];
f_uxx_n25 = [9.25e5,1.58e3,1.80e4,3.06];
f_uxx_n50 = [1.23e7,2.90e4,5.86e4,4.27e1];
f_uxx = [f_uxx_n0;f_uxx_n1;f_uxx_n5;f_uxx_n10;f_uxx_n25;f_uxx_n50];

plot_RMSE(f_u,f_ut,f_ux,f_uxx,'\textbf{Fishers}');

%% Extension 1: comparing different ANN structures for advection diffusion

% 1000 neurons (original)

u_1000 = u(:,4);
ut_1000 = ut(:,4);
ux_1000 = ux(:,4);
uxx_1000 = uxx(:,4);

% 500 neurons

u_500 = [6.53e-4,3.83e-4,1.33e-3,1.05e-3,6.29e-3,9.14e-2];
ut_500 = [4.01e-1,4.67e-1,2.44,2.66,3.89,1.01e1];
ux_500 = [1.06e-2,2.61e-2,7.08e-2,2.93e-2,9.17e-2,4e-1];
uxx_500 = [7.25e-1,4.64e-1,1.01,9.44e-1,6.50e-1,3.23];

% 2000 neurons

u_2000 = [3.57e-4,4.13e-4,7.94e-4,6.46e-4,4.24e-3,7.26e-2];
ut_2000 = [2.54e-1,9.58e-2,1.88e-1,6.79e-1,1.26,2.88];
ux_2000 = [1.22e-2,7.04e-3,3.69e-2,4.40e-2,4.28e-2,1.46e-1];
uxx_2000 = [4.17e-1,6.34e-1,6.65e-1,9.14e-1,7.76e-1,1.44];

u_3000 = [4.45e-4,3.76e-4,3.75e-4,9.79e-4,6.93e-3,8.48e-2];
ut_3000 = [3.40e-1,5.79e-1,1.43e-1,7.03e-1,1.48,3.13];
ux_3000 = [2.09e-2,3.40e-2,1.08e-2,2.35e-2,4.89e-2,1.28e-1];
uxx_3000 = [4.41e-1,5.12e-1,5.80e-1,8.76e-1,5.87e-1,1.63];

neur_u = [u_500;u_1000';u_2000;u_3000]';
neur_ut = [ut_500;ut_1000';ut_2000;ut_3000]';
neur_ux = [ux_500;ux_1000';ux_2000;ux_3000]';
neur_uxx = [uxx_500;uxx_1000';uxx_2000;uxx_3000]';

plot_normRMSE_heatmap(neur_u,neur_ut,neur_ux,neur_uxx,'neurons','Advection Diffusion');

%% Extension 3: comparing ANN data denoising results for reduced data

u_0_ANN = f_u(:,4)';
u_10_ANN = [3.12e-3,3.16e-3,2.95e-3,3.67e-3,8.53e-3,5.55e-2];
u_50_ANN = [1.171e-3,1.76e-3,1.66e-3,2.22e-3,9.79e-3,7.18e-2];
u_80_ANN = [5.82e-2,1.25e-3,4.55e-2,2.24e-2,4.09e-2,1.81e-1];
u_vec_ANN = [u_0_ANN;u_10_ANN;u_50_ANN;u_80_ANN]';

ut_0_ANN = f_ut(:,4)';
ut_10_ANN = [3.05e-2,4.17e-2,1.11e-1,7.12e-2,3.38e-2,8.28e-1];
ut_50_ANN = [1.70e-1,2.49e-1,2.21e-1,9.42e-2,4.37e-1,3.69e-1];
ut_80_ANN = [4.49e-1,2.24e-1,3.16e-1,1.52e-1,3.68,3.89];
ut_vec_ANN = [ut_0_ANN;ut_10_ANN;ut_50_ANN;ut_80_ANN]';

ux_0_ANN = f_ux(:,4)';
ux_10_ANN = [2.49e-2,3.67e-2,7.31e-2,8.74e-2,2.46e-2,1.13];
ux_50_ANN = [2.46e-2,2.08e-2,2.23e-2,1.52e-1,3.98e-1,4];
ux_80_ANN = [6.37e-1,7.24e-2,3.41e-1,3.74e-1,6.14,3.11];
ux_vec_ANN = [ux_0_ANN;ux_10_ANN;ux_50_ANN;ux_80_ANN]';

uxx_0_ANN = f_uxx(:,4)';
uxx_10_ANN = [3.06e-1,3.55e-2,1.11e-1,1.29e-1,4.10e-1,2.66e1];
uxx_50_ANN = [1.44e-1,2.19e-1,5.73e-1,9.59e-1,1.71e1,1.54e2];
uxx_80_ANN = [3.37e-1,2.96e-2,4.23e-1,6.91e-1,2.3,3.64];
uxx_vec_ANN = [uxx_0_ANN;uxx_10_ANN;uxx_50_ANN;uxx_80_ANN]';

plot_normRMSE_heatmap(u_vec_ANN,ut_vec_ANN,ux_vec_ANN,uxx_vec_ANN,'reduction','Fishers');

%% PDEFIND Fishers

FD = [0.001764,0.3375212,1,1,1,1];
LCVSP = [0.00845,0.06305,0.3623988,1,1,1];
LNCVSP = [0.0075555,0.0143247,0.338905,0.341047133,1,1];
ANN = [0.0236504,0.055789,0.023968,0.044709,0.0573424,0.69];
fisher_score = [FD;LCVSP;LNCVSP;ANN]';

plot_PDEFIND_bar(fisher_score,'\textbf{Fishers}');

