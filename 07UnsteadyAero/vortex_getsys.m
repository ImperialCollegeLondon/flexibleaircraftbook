% Function vortex_getsys
%  Compute state-space matrices of the vortex solution on the airfoil.
%  Inputs: 
%    Nb: Number of panels on the airfoil
%    Nw: Number of panels on the wake
%    xa: Relative coordinates from the vortex locations to the reference point
%        for moment calculations.
% 
function [A,B,C,D]=vortex_getsys(Nb,Nw,xa)

% Define the value of the dissipation factor. These numbers have been
% obtained to enforce that the largest eigenvalue is 0.001, but it makes
% no actual difference to most results. It can be set to zero for bode
% plots.
epsilon=2.5e-4;
if Nw==Nb*30
    if Nb==100
        epsilon=0.2074e-4;
    elseif Nb==20
        epsilon=1.0193e-4;
    elseif Nb==8
        epsilon=2.5439e-4;
    end
end
if Nb==100
    if Nw==Nb*5
        epsilon=0.2224e-4;
    elseif Nw==Nb*15
        epsilon=0.2100e-4;
    elseif Nw==Nb*30
        epsilon=0.2074e-4;
    elseif Nw==Nb*100
        epsilon=0.2122e-4;
    elseif Nw==Nb*200
        epsilon=0.2426e-4;
    else
        epsilon=0.2e-4;
    end
end

%% Equations are written as: E x_a(n+1) = F x_a(n) + G w(n+1)
E=sparse(Nb+Nw,Nb+Nw);
F=sparse(Nb+Nw,Nb+Nw);
G=sparse(Nb+Nw,Nb);

% Non-penetration boundary conditions. 
for j=1:Nb
    for k=1:Nb+Nw
       % E(j,k)=1/(2*pi*(x(k)-xi(j)));
        E(j,k)=Nb/(2*pi*(k-j-1/2));
    end
end
G(1:Nb,1:Nb)=-eye(Nb);

% Kutta condition.
E(Nb+1,1:Nb+1)=1;
F(Nb+1,1:Nb)  =1;

% Wake convection.
for j=2:Nw
  E(Nb+j,Nb+j)  =1;
  F(Nb+j,Nb+j-1)=1;
end
F(Nb+j,Nb+j)=1-epsilon;  % Accummulate/dissipate vorticity in last panel.

% State matrix.
A=full(E\F);

% Input matrix.
B=full(E\G);

%% Output equations
C=zeros(2,Nb+Nw);
% Lift:
for j=1:Nb
  C(1,j)=2-2*(Nb-j+1);
  Y(1,j)=2*(Nb-j+1);
end

% Moment:
for j=1:Nb
   beta=sum(xa(j:Nb));
   C(2,j)=-2*xa(j)+2*beta;
   Y(2,j)=-2*beta;
end

C=C+Y*A(1:Nb,:);
D=Y*B(1:Nb,:);
% eof