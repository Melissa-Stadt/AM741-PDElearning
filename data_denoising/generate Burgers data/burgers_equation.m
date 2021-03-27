%u_t = eta u_xx + uu_x

clear all; close all; 

% I am going to use pdepe
x = linspace(-1,1,201);
t = linspace(0,1,101);
dx = x(2)-x(1);
dt = t(2) - t(1);
% initialize
eta = 0.6;
m=0;
eqn=@(x,t,u,dudx) burgersPDE(x,t,u,dudx, eta);
ic=@(x) burgerIC(x);
sol=pdepe(m,eqn,ic,@burgerBC,x,t);

U = sol(:,:,1);
% initialize
U_t = zeros(length(t)-2,length(x)-2);
U_x = U_t;
U_xx = U_x;

% at each time step take the U_x values
for iter =2:(size(U,1)-1)
    % x values at time t(iter)
    U_it = U(iter,:);
    [U_x_it, U_xx_it] = FiniteDifferences_x(U_it,dx);
    U_x(iter-1,:) = U_x_it;
    U_xx(iter-1,:) = U_xx_it;
end %for

U_t = eta.*U_xx + U(2:end-1, 2:end-1).*U_x;
% at each spatial step take the U_t values
% for x_iter = 2:(size(U,2)-1)
%     % time values at each spatial value
%     U_time = U(:,x_iter);
%     U_t_it = FiniteDifference_t(U_time,dt);
%     U_t(:,iter-1) = U_t_it;
% end %for

Umax = max(max(U));
Umin = min(min(U));


% save data
filename = 'burgers_groundtruth.mat';
save(filename, 'U','U_t', 'U_x', 'U_xx', 't', 'x')


n_it = length(t)
x_ends = x(2:end-1);
figure(1)
hold on
for iter =1:n_it-2
    disp(iter)
    Uit = U(iter,:);
    U_xit = U_x(iter,:);
    U_tit = U_t(iter,:);
    U_xxit = U_xx(iter,:);
    % plot results
    subplot(2,2,1)
    plot(x,Uit)
    xlabel('x')
    ylabel('U')
    ylim([Umin, Umax])
    subplot(2,2,2)
    plot(x_ends,U_xit)
    xlabel('x')
    ylabel('U_x')
    subplot(2,2,3)
    plot(x_ends,U_tit)
    xlabel('x')
    ylabel('U_t')
    subplot(2,2,4)
    plot(x_ends,U_xxit)
    xlabel('x')
    ylabel('U_{xx}')
    pause(0.0001)
end
hold off

figure(10)
surf(x,t,U)
title('burgers')
xlabel('x')
ylabel('t')
zlabel('solution u')

figure(99)
plot(x, burgerIC(x))
title('initial condition')

[uout,duoutdt] = pdeval(0,t,sol(:,2,1),0.5)

    
%-----------------------------------------------------------------
% functions used
%----------------------------------------------------------------

%PDE function
function [c,f,s] = burgersPDE(x,t,u,dudx, eta)
c = 1;
f=eta*dudx;
s=u*dudx;
end %PDE

% initial condition
function u0 = burgerIC(x)
a=1.0;
b=0;
c=0.5;
u0 = a/(c*sqrt(2*pi))*exp(-0.5*(x-b).^2/c.^2);
end %IC

% boundary condtions
function [pl, ql, pr, qr]=burgerBC(xl,ul,xr,ur,t)
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

% dt is the timestep
% U_it is the time data
function U_t_it = FiniteDifference_t(U_time, dt)
    n = length(U_time);
    U_t_it = zeros(n-2,1);
    for ii = 2:n-1
        jj = ii-1; % index for derivative vectors
        U_t_it(jj) = (U_time(ii+1)-U_time(ii-1))/(2*dt);
    end %for
end %function