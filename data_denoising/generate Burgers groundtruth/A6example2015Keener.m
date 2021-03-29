% This is an example from A.6 of the 2015KeenerPDEsofBiology_coursenotes
% A.6 Diffusion equation with growth/decay solution via Crank-Nicolson

clear all; close all

L=1000; % length of the domain
d = 0.01; %diffusion coefficient

N=50; %number of spatial grid points
M=130; %number of time steps
h = L/N;
dt = h/2;
K=-1; % pick K<0 for decay, K>0 for growth

if K>0
    kp=K;
    km=0;
else
    kp = 0;
    km=K;
end % if

sc=100;
scal=d*dt/h^2;

X = h*(1:N)';
V =zeros(N,1);

V(N/2) = 1; %iniitialize V

kstep = 10;
t=0;

% uses Crank Nicolson to solve the diffusion equation
% uses explicit Euler for decay, implicit Euler for growth
% set up matrix

Atm=(1 + scal)*ones(N,1)-kp;
Atm(1,1) = 1 + scal/2-kp;
Atm(N,1) = 1+scal/2-kp;

Am=diag(Atm)-diag(scal/2*ones(N-1,1),1)-diag(scal/2*ones(N-1,1),-1);

% algorithm

for n = 2:M
    for k = 1:kstep
        F=vel*(V-[V(1);V(1:N-1)]);
    end % k for
end % n for