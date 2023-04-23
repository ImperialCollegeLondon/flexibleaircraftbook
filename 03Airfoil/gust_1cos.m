%
% gust_1cos: Compute 1-cos response for 1-DoF airfoil, using Jones's
%            approximation for Theodorsen's lift deficiency function.
%
% Dependencies:
%    sears.m: Analytical expression for Sears's function.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all 

% Coefficients in Jones's approximation for C(ik)
a_1=0.165;
a_2=0.335;
b_1=0.0455;
b_2=0.3;


% Obtain RFA approximation to Sears
k=0:0.02:50;
sysse=frd(sears(k),k);
ssys5=fitmagfrd(sysse,5);

mu=5;

for H=5:5:50;  % Gradient in semichords.
    
  % Create 1-cos gust velocity profile.
  t=0:H/100:10*H;
  wg(1:200)=0.5*(1-cos(pi*t(1:200)/H));
  wg(201:1001)=0;

  subplot(3,1,1)
   plot(t,wg,'k','LineWidth',2), hold on
   ylabel('w_g/V_\infty','FontSize',14,'FontWeight','bold')
   axis([0 120 0 1]) 
  
  % Obtain lift due to gust velocity.
  [CLg, tg]=lsim(ssys5,wg,t)  % Lift divided by 2\pi.

  subplot(3,1,2)
   plot(tg,CLg,'k','LineWidth',2), hold on
   ylabel('C_L^g/(2\pi)','FontSize',14,'FontWeight','bold')
   axis([0 120 0 1])

  % Equations of motion for the 1-dof aerofoil. 
  A_1=[[4*mu+1 2*a_1 2*a_2];[0 1 0];[0 0 1]];
  A_2=[-2 0 0; -1 -b_1 0; -1 0 -b_2];
  A=inv(A_1)*A_2;
  B0=[1; 0; 0];
  B=inv(A_1)*B0;
  
  % Output is the acceleration.
  C=A(1,:);
  D=B(1,1);
  sysgust=ss(A,B,C,D);

  [y, t]=lsim(sysgust,CLg,tg)
  subplot(3,1,3)
    plot(t,y,'k','LineWidth',2)
    hold on
    axis([0 120 -0.02 0.04])
    xlabel('Reduced time, s','FontSize',14,'FontWeight','bold')
    ylabel('d\nu/ds','FontSize',14,'FontWeight','bold')
end

