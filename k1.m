clear all

N=10;

L=1.5;
H=0.5;

T1= @(x) 5*(x/L-1) + 15*cos(pi*x/L);
T2=10;
T3=15;
% dT/dx=0 at boundary 4 => a_W=0

k1=20;
k2=0.01;
b=0; % source

% bd 1: y=0
% bd 2: x=L
% bd 3: y=H
% bd 4: x=0

xface=linspace(0, L, N-1);
yface=linspace(0, H, N-1);

xnode=zeros(N, 1);
ynode=zeros(N, 1);
