% gust_1cos.m
% 
%         Compute 1-cos response for 1-DoF airfoil. Section 3.6.2 of
%         Palacios & Cesnik's book.
%
% Dependencies:
%    theodorsen_ajbj.m: aj & bj coefficients of RFA to Theodorsen.
%    sears.m: Analytical expression for Sears's function.
%
% Copyright, Rafael Palacios, May 2023
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all 

% RFA of Theodorsen's Lift deficiency function.
[a,b]=theodorsen_ajbj(2);
a_1=a(1); a_2=a(2); b_1=b(1);b_2=b(2);

% Uncomment to use Jones's coefficients instead.
% a_1=0.165; a_2=0.335; b_1=0.0455; b_2=0.3;

% Obtain RFA approximation to Sears of order 5.
ssys5=sears_rfa(5);

mu=5;     % mass parameter.
dH=100;   % Discretization of the gust profile.
NumH=10;  % Number of gusts lengths that the simulation runs for.

for H=5:5:50;  % Gradient in semichords.
    
  % Create 1-cos gust velocity profile.
  t=0:H/dH:NumH*H;
  wg(1:2*dH)=0.5*(1-cos(pi*t(1:2*dH)/H));
  wg(2*dH+1:NumH*dH+1)=0;

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

