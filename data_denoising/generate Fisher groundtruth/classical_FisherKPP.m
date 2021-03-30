%--------------------------------------------------------
% This file generates the new_fisher_groundtruth.mat
% data for given classical Fisher-KPP equation data
% Will compute the solution and the U_t, U_xx, and U_x
% derivative values using Finite Differences

% Options to simulate the output and compare to the 
% original groundtruth data file are there
%-------------------------------------------------------

clear all; close all;

%------------------------------------------------------------
% Begin user input
%------------------------------------------------------------
% cell or original gaussian
% if want cell then set to 1, else gaussian profile
cell = 0;

% if want to show simulation of FisherKPP set show_sim to 1
% NOTE: this is long for cell since time goes to 200
show_sim = 1;

% if want to show initial condition
show_IC = 0;

% test compares to the original grountruth data
% set to 1 if want to compare
% additional flag options in script for simulations
% default shows all
test = 0;
%------------------------
% if testing, load original data
%------------------------
if test == 1
    % old file for testing
    load('fisher_groundtruth.mat')
    Dold = D;
    Kold = K;
    Uold = U;
    U_told = U_t;
    U_xold = U_x;
    U_xxold = U_xx;
    r_old = r;
    t_old = t;
    x_old = x;
    flag = 1; % full time simulation, note time series must be same
    flag1 = 1; % show surface plots of original data and new
    flag2 = 1; % show initial conditions
else
    flag = 0;
    flag1 = 0;
    flag2 = 0;
end %testing flags


%----------------------------------------
%Parameters for Fisher KPP
%---------------------------------------
% cell parameters for initial cell density 10,000 from BINNs paper
if cell == 1
    D = 309.7;
    r = 0.0437;
    K = 1.7e03;
else
    D = 0.06; % original D = 0.02
    K = 1; % original K = 1
    r = 8; % original r = 10
end % if cell
params = {D, K, r};


% mesh
if cell == 1
    x = linspace(0, 2, 201);
    t_end = 200;
    t_vals = 20*t_end;
    t = linspace(0, t_end, t_vals);
else
    x = linspace(-1,1,201);
    t = linspace(0,1,101);
end

dx = x(2)-x(1);
dt = t(2) - t(1);
%-------------------------------------
% end user input
%-------------------------------------

% solve PDE
m=0;
eqn=@(x,t,u,dudx) FisherPDE(x,t,u,dudx, params);
if cell == 1
    ic=@(x) cell_densityIC(x); % cell density profile
else
    ic=@(x) FisherIC(x); % this is for the original model
end

sol=pdepe(m,eqn,ic,@FisherBC,x,t);

U = sol(:,:,1);

% at each time step take the U_x values
for iter =2:(size(U,1)-1)
    % x values at time t(iter)
    U_it = U(iter,:);
    [U_x_it, U_xx_it] = FiniteDifferences_x(U_it,dx);
    U_x(iter-1,:) = U_x_it;
    U_xx(iter-1,:) = U_xx_it;
end %for

%at each spatial step take the U_t values
for x_iter = 2:(size(U,2)-1)
    % time values at each spatial value
    U_time = U(:,x_iter);
    U_t_it = FiniteDifference_t(U_time,dt);
    U_t(:,x_iter-1) = U_t_it;
end %for

%-------------------------------------
% save data
%-------------------------------------
if cell == 1
    filename = 'cell_fisherKPP_groundtruth.mat';
else
    filename = 'new_fisher2_groundtruth.mat';
end
save(filename, 'U', 'U_t', 'U_x', 'U_xx', 't', 'x', 'K', 'D', 'r')

if show_sim == 1
    n_it = length(t)
    x_ends = x(2:end-1);
    figure(1)
    hold on
    for iter =1:n_it-2
        disp(t(iter))

        % values at iteration
        Uit = U(iter,:);
        U_xit = U_x(iter,:);
        U_xxit = U_xx(iter,:);
        U_tit = U_t(iter,:);

        % plot results
        subplot(2,2,1)
        % U
        plot(x,Uit, 'b-', 'linewidth', 2)
        hold on
        xlabel('x')
        ylabel('U')
        hold off


        subplot(2,2,2)
        % U_x
        plot(x_ends,U_xit, 'b-', 'linewidth', 2)
        hold on
        xlabel('x')
        ylabel('U_x')
        hold off


        subplot(2,2,3)
        plot(x_ends,U_tit, 'b-','linewidth',2)
        hold on
        xlabel('x')
        ylabel('U_t')
        hold off

        subplot(2,2,4)
        plot(x_ends,U_xxit, 'b-', 'linewidth', 2)
        hold on
        xlabel('x')
        ylabel('U_{xx}')
        hold off
        pause(0.0001)
    end
    hold off
end %show_sim

if show_IC == 1
   figure(50)
   plot(x, U(1,:))
   title('initial condition')
   xlabel('x')
   ylabel('U')
end

%----------------------------------------
% below are figure options for testing
%----------------------------------------

if flag2 == 1
    figure(20)
    surf(x,t,U)
    title('Fisher KPP')
    xlabel('x')
    ylabel('t')
    zlabel('U new solution')

    figure(21)
    surf(x_old,t_old,Uold)
    title('original fisher KPP')
    xlabel('x')
    ylabel('t')
    zlabel('U original')
end %flag2



if flag1 == 1
    figure(22)
    % initial condition
    plot(x_old, Uold(1,:), 'linewidth', 2)
    hold on
    plot(x, FisherIC(x), 'linewidth', 2)
    xlabel('x')
    ylabel('U_0')
    legend('original U0', 'FisherIC')
    hold off
end % flag1

if flag2 == 1
    figure(23)
    surf(x,t,U)
    title('Fisher KPP')
    xlabel('x')
    ylabel('t')
    zlabel('U new solution')

    figure(24)
    surf(x_old,t_old,Uold)
    title('original fisher KPP')
    xlabel('x')
    ylabel('t')
    zlabel('U original')
end %flag2

if flag == 1
    n_it = length(t)
    x_ends = x(2:end-1);
    figure(99)
    hold on
    for iter =1:n_it-2
        disp(iter)
        % original values
        Uit = Uold(iter,:);
        U_xit = U_xold(iter,:);
        U_tit = U_told(iter,:);
        U_xxit = U_xxold(iter,:);
        
        % new values
        Unewit = U(iter,:);
        U_xnewit = U_x(iter,:);
        U_xxnewit = U_xx(iter,:);
        U_tnewit = U_t(iter,:);
        
        % plot results
        subplot(2,2,1)
        % U
        plot(x,Uit, 'b-', 'linewidth', 2)
        hold on
        plot(x,Unewit, 'r--', 'linewidth', 2)
        xlabel('x')
        ylabel('U')
        %legend('U original', 'U_{new}')
        hold off
        
        
        subplot(2,2,2)
        % U_x
        plot(x_ends,U_xit, 'b-', 'linewidth', 2)
        hold on
        plot(x_ends,U_xnewit, 'r--', 'linewidth', 2)
        xlabel('x')
        ylabel('U_x')
        hold off
        
        
        subplot(2,2,3)
        plot(x_ends,U_tit, 'b-','linewidth',2)
        hold on
        plot(x_ends, U_tnewit, 'r--', 'linewidth', 2)
        xlabel('x')
        ylabel('U_t')
        hold off
        
        subplot(2,2,4)
        plot(x_ends,U_xxit, 'b-', 'linewidth', 2)
        hold on
        plot(x_ends, U_xxnewit, 'r--', 'linewidth', 2)
        xlabel('x')
        ylabel('U_{xx}')
        hold off
        pause(0.001)
    end
    hold off
end %flag if



%-------------------------------------------------------------
% functions used in file
%-------------------------------------------------------------

% PDE function: classical Fisher-KPP
function [c,f,s]=FisherPDE(x,t,u,dudx,params)
D = params{1};
K = params{2};
r = params{3};
c=1;
f=D*dudx;
s =r*u*(1-u/K);
end % PDE

% initial condition
function u0 = FisherIC(x)
a = 0.0075;
b = 0;
c = 0.03;
u0 = a/(c*sqrt(2*pi))*exp(-0.5*(x-b).^2/c.^2);
end %IC

%boundary conditions
function [pl, ql, pr, qr]=FisherBC(xl,ul,xr,ur,t)
pl =ul;
ql=0;
pr=ur;
qr=0;
end %BC

% finite differences for computing derivative
% computes the first and second derivative of data using 
% the finite difference method
% Inputs:
% dx is the spatial step size
% U is the spatial data
function [U_x_it, U_xx_it] = FiniteDifferences_x(U_it, dx)
    n = length(U_it);
    U_x_it = zeros(n-2, 1);  % will not have the endpoints
    U_xx_it = zeros(n-2, 1);
    for ii = 2:n-1
        jj = ii - 1; % index for the derivative vectors
        U_x_it(jj) = (U_it(ii+1)-U_it(ii-1))/(2*dx);
        U_xx_it(jj) = (U_it(ii+1)-2*U_it(ii) + U_it(ii-1))/dx.^2;
    end % for
end % FiniteDifference

function U_t_it = FiniteDifference_t(U_time, dt)
    n = length(U_time);
    U_t_it = zeros(n-2,1);
    for ii = 2:n-1
        jj = ii-1; % index for derivative vectors
        U_t_it(jj) = (U_time(ii+1)-U_time(ii-1))/(2*dt);
    end %for
end %function

