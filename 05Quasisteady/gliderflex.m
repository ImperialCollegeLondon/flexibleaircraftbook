%
% gliderflex.m: Stability of the simple glider with a flexible fuselage.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Define constants in the problem.
m=318;
Iyy=432;
g=9.81;
St=1.14;
S=6.11;
b=15;
lw=0.3;
lt=4.63;
lwlt=lw+lt;
lb=lwlt/2;

% Aerodynamic coefficients.
rho=1.112; % 1000 m
e=0.8;
CLaw=5.55;
CLat=4.30;
CLdt=0.45;

% Fligth conditions.
Vinf=30;
theta0=-2*pi/180;
ref=1/2*rho*(Vinf^2)*St*lt;

% Loop through the stiffness of the fuselage spring.
k=0;
for kb=[1.87 2:2:20 50 10000]*ref
    k=k+1;
    lt=lwlt-lw;
    qinf=1/2*rho*Vinf^2;
    CLa=CLaw+(St/S)*CLat;
    A_r=[2*g*sin(theta0)/Vinf    g*cos(theta0)/Vinf*(1-(2*S*CLa)/(pi*b^2*e))  0                              -g*cos(theta0);
        -2*g*cos(theta0)/Vinf    g*sin(theta0)/Vinf-rho*Vinf*S*CLa/(2*m)      Vinf-rho*Vinf*St*lt*CLat/(2*m) -g*sin(theta0);
         0                        rho*Vinf*(S*lw*CLaw-St*lt*CLat)/(2*Iyy)      -rho*Vinf*St*lt^2*CLat/(2*Iyy) 0;
         0                        0                                            1                              0];

    A_rs=qinf*St*CLat/m*[0; 1; lt*m/Iyy; 0];
    A_sr=1/2*rho*Vinf*CLat*St*lb*[0 1 lt 0];
    A_s =-qinf*St*lb*CLat;
    
    A=A_r+A_rs*A_sr*(1/(kb-A_s))

    % Obtain the eigenvalues of the aircraft dynamics.
    [V, D]=eig(A);
    s_long(:,k)=diag(D);
    
end

% Plot the eigenvalues.
  plot(real(s_long(:,1:end-1)),imag(s_long(:,1:end-1)),'kx','MarkerSize',8)
  hold on
  plot(real(s_long(:,end)),imag(s_long(:,end)),'ko','MarkerSize',8)
  
  axis([-3.5 0 0 4])
  xlabel('Re(\lambda)','FontSize',14)
  ylabel('Im(\lambda)','FontSize',14)
  grid on
% eof
