%u_t = eta u_xx + uu_x

clear all; close all; 

% I am going to use pdepe
x = linspace(-1,1,200);
t = linspace(0,1,100);
% initialize
eta = 1;
m=0;
eqn=@(x,t,u,dudx) burgersPDE(x,t,u,dudx, eta);
ic=@(x) burgerIC(x,eta);
sol=pdepe(m,eqn,ic,@burgerBC,x,t);

u = sol(:,:,1);
filename = 'burgers_groundtruth.mat';
save(filename, 'u','t', 'x')
% figure(1)
% surf(x,t,u)
% title('burgers')
% xlabel('x')
% ylabel('t')
% zlabel('solution u')
% figure(1)
% plot(x, burgerIC(x,1))
%------------------
% functions used
%-------------------

%PDE function
function [c,f,s] = burgersPDE(x,t,u,dudx, eta)
c = 1;
f=eta*dudx;
s=u*dudx;
end %PDE

% initial condition
% Gaussian profile
function u0 = burgerIC(x, eta)
a=0.1;
b=0;
c=2;
u0=a*exp(-(x-b).^2/(2*c.^2));
end %IC

% boundary condtions
function [pl, ql, pr, qr]=burgerBC(xl,ul,xr,ur,t)
pl =ul;
ql=0;
pr=ur;
qr=0;
end %BC