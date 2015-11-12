% Group 28 - Sebastian Ståhl and Kohei Obara

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

for i=2:N-1
    xnode(i)=xface(i-1) + (xface(i)-xface(i-1))/2;
    ynode(i)=yface(i-1) + (yface(i)-yface(i-1))/2;
end

xnode(end)=L;
ynode(end)=H;

% k=k2=0.01 in 0.7<x<1.1, 0.3<y<0.4
% k=k1=20   otherwise

% k contains info of k in nodes
k=k1*ones(N, N);

xlowlim=find((xnode-0.7)>0, 1);
xuplim=find((xnode-1.1)>0, 1)-1;

ylowlim=find((ynode-0.3)>0, 1);
yuplim=find((ynode-0.4)>0, 1)-1;

k(ylowlim:yuplim, xlowlim:xuplim)=k2;

T=zeros(N,N);

% enforce b.c.'s
T(:, end)=T2;
T(end, :)=T3;
T(1, :)=T1(xnode);

iter=0;


while(iter<1)
    
    R=0; % residual
    
    
    for i=2:N-1
        
        Delta_x=xface(i)-xface(i-1);
        dx_e=xnode(i+1)-xnode(i);
        dx_w=xnode(i)-xnode(i-1);
        
        for j=2:N-1
            
            Delta_y=yface(j)-yface(j-1);
            dx_n=ynode(j+1)-ynode(j);
            dx_s=ynode(j)-ynode(j-1);
            
            % calculate k at faces through interpolation
            fx=0.5*Delta_x/dx_e;
            k_e=fx*k(i+1,j) + (1-fx)*k(i,j);
            
            fx=0.5*Delta_x/dx_w;
            k_w=fx*k(i-1,j) + (1-fx)*k(i,j);
            
            fy=0.5*Delta_y/dx_s;
            k_s=fy*k(i,j-1) + (1-fy)*k(i,j);
            
            fy=0.5*Delta_y/dx_n;
            k_n=fy*k(i,j+1) + (1-fy)*k(i,j);
            
            % calculate contribution coeffs
            a_E=k_e*Delta_y/dx_e;
            a_W=k_w*Delta_y/dx_w;
            a_S=k_s*Delta_x/dx_s;
            a_N=k_n*Delta_x/dx_n;
            a_P=a_E+a_W+a_S+a_N;
            
            % dT/dx=0 at boundary 4 => a_W=0
            if(i==1)
                a_W=0;
            end
            
            T_E=T(j, i+1);
            T_W=T(j, i-1);
            T_S=T(j-1, i);
            T_N=T(j+1, i);
            
            T(j, i)=(a_E*T_E+a_W*T_W+a_S*T_S+a_N*T_N)/a_P;
            
            R=R+abs(a_E*T_E+a_W*T_W+a_S*T_S+a_N*T_N - T(j,i)*a_P);
            
        end
    end
    
    iter=iter+1;
end

%
figure(1)
clf
hold on
surf(xnode, ynode, T)
xlabel('x')
ylabel('y')
view([-7 26])
