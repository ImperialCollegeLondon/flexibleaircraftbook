%
% glider.m: Stability of a simple glider.
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

% Aerodynamic constants.
rho=1.112; % Density for 1000 m
e=0.8;     % Oswall efficiency factor.
CLaw=5.55;
CLat=4.30;
CLdt=0.45;

% Flight conditions.
Vinf=30;
theta0=-2*pi/180;
alphat0=-5*pi/180;

%% Trim the aircraft.
% load factor
n_trim=cos(theta0);
A=[[S*CLaw+St*CLat St*CLdt]; ...
   [lw*S*CLaw-lt*St*CLat -lt*St*CLdt]];
B=[n_trim*m*g/(0.5*rho*Vinf^2)-St*CLat*alphat0; ...
   lt*St*CLat*alphat0];

X=inv(A)*B;
alpha_trim=X(1)*180/pi
delta_trim=X(2)*180/pi


% Choose a variable (lw here) and loop through it.
k=0;
ref=lw;
for lw=[ref 0.8*ref:ref/10:ref-ref/10 ref+ref/10:ref/10:1.2*ref]
    k=k+1;
    lt=lwlt-lw;
    qinf=1/2*rho*Vinf^2;
    CLa=CLaw+(St/S)*CLat;
    A=[2*g*sin(theta0)/Vinf    g*cos(theta0)/Vinf*(1-(2*S*CLa)/(pi*b^2*e))  0                              -g*cos(theta0);
      -2*g*cos(theta0)/Vinf    g*sin(theta0)/Vinf-rho*Vinf*S*CLa/(2*m)      Vinf-rho*Vinf*St*lt*CLat/(2*m) -g*sin(theta0);
       0                        rho*Vinf*(S*lw*CLaw-St*lt*CLat)/(2*Iyy)      -rho*Vinf*St*lt^2*CLat/(2*Iyy) 0;
       0                        0                                            1                              0];

    B=-[0; qinf*St*CLdt/m; qinf*St*lt*CLdt/Iyy; 0];

    [V, D]=eig(A);
    s_long(:,k)=diag(D);
    
end

% Plot the eigenvalues as a funciton of the variable.
  plot(real(s_long(:,1)),imag(s_long(:,1)),'ro')
  hold on
  plot(real(s_long(:,2:end)),imag(s_long(:,2:end)),'kx')
  axis([-3.5 0 0 4.5])
  xlabel('Re(\lambda)','FontSize',14)
  ylabel('Im(\lambda)','FontSize',14)
  grid on
% eof
