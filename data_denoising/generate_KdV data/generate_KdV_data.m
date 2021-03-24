% generate two soliton KdV equation

clear all; close all; clc

% initialize

% spatial range
x_start = -1;
x_end = 1;
% spatial step
dx = 0.01;
% spatial grid
% periodic conditions u(-1) = u(1)
x = (x_start:dx:x_end)';
N = length(x);

% temporal range
t_start = 0;
t_end = 1;
% time step
% Note: Goda is unconditionally stable
dt = 0.01;
%temporal grid
t = (0:dt:t_end);
n_it = length(t);

% initialize
U = zeros(length(t),length(x));
U_t = zeros(length(t)-2,length(x)-2);
U_x = U_t;
U_xxx = U_x;

%---------------------------
% initial condition
%--------------------------
% two soliton initial condition
v_old = -(8./((cosh(2.*x+8).^2))) - (2./((cosh(x)).^2));
%---------------------------
% main iteration loop
%---------------------------
for iter = 1:n_it
    v_new = Goda(v_old, dt, dx);
    % plot output
%     plot(x,v_new, '*b-')
%     xlabel('x')
%     ylabel('v')
%     xlim([x_start, x_end])
%     ylim([-2.5, -0.5])
%     pause(0.00001)
    U(iter, :) = v_new;
    % set up for next time step
    v_old = v_new;
end % end main iter loop

for iter = 2:(size(U,1) - 1)
    U_it = U(iter, :);
    [U_x_it, U_xxx_it] = FiniteDifferences_x(U_it, dx);
    U_x(iter-1, :) = U_x_it;
    U_xxx(iter-1,:) = U_xxx_it;
end % for

for x_iter = 2:(size(U,2) - 1)
    U_it = U(:, x_iter);
    U_t_it = FiniteDifference_t(U_it,dt);
    U_t(:, iter-1) = U_t_it;
end % for
% filename to save data in
filename = 'KdVgroundtruth.mat';
% format correctly for data files
x = x'; 
t = t';
save(filename, 'U', 'U_t', 'U_x', 'U_xxx', 't', 'x')


%-----------------------------------------------------------
% functions 
%-----------------------------------------------------------
% finite differences for computing derivative
% computes the p-th derivative of data using 
% the finite difference method
% Inputs:
% p <= 3 is the derivative to take
% dx is the spatial step size
% U is the spatial data
function [U_x_it, U_xxx_it] = FiniteDifferences_x(U_it, dx)
    n = length(U_it);
    U_x_it = zeros(n-2, 1);  % will not have the endpoints
    U_xxx_it = zeros(n-2, 1);
    for ii = 2:n-1
        jj = ii - 1; % index for the derivative vectors
        U_x_it(jj) = (U_it(ii+1)-U_it(ii-1))/(2*dx);
        if ii == 2
            U_xxx_it(jj) = (-2.5*U_it(1) + 9*U_it(2) - 12*U_it(3) + 7*U_it(4) - 1.5*U_it(5))/dx.^3;
        elseif ii == n-1
            U_xxx_it(jj) = (2.5*U_it(n) - 9*U_it(n-1) + 12*U_it(n-2) - 7*U_it(n-3) + 1.5*U_it(n-4))/dx.^3;
        else
            U_xxx_it(jj) = (U_it(ii+2)/2-U_it(ii+1)+U_it(ii-1)-U_it(ii-2)/2)/dx.^3;
        end
    end % for
end % FiniteDifference

% dt is the timestep
% U_it is the time data
function U_t_it = FiniteDifference_t(U_it, dt)
    n = length(U_it);
    U_t_it = zeros(n-2,1);
    for ii = 2:n-1
        jj = ii-1;
        U_t_it(jj) = (U_it(ii+1)-U_it(ii-1))/(2*dt);
    end %for
end %function








