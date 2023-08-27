% vortex_bode.m
%
%  Compare transfer functions between alpha/beta and CL/CM using a 
%  state-space vortex description and the analytical formulas.
%
% Dependencies:
%    theodorsen.m: Analytical expression for Theodorsen's lift deficiency
%                  function.
%
% It supports Section 7.2.4 in Palacios & Cesnik (CUP, 2023)
%    https://www.cambridge.org/9781108420600
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: May 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all, close all
addpath('../03Airfoil/')

%% Define problem parameters. All dimensions correspond to unit chord.
Nb=20;                         % Number of segments along the aerofoil
Nw=Nb*30;                      % Number of chordlengths in wake
dx=1/Nb;                       % Nondimensional panel length
x=dx/4:dx:1+Nw/Nb-((3*dx)/4);  % Coordinates of vortices (aerofoil/wake).
xi=3*dx/4:dx:1-dx/4;           % Coordinates of collocation points

x0=1/4;                        % Pitch about the quarter chord.

wbode=0.001:0.001:2;           % Frequencies for the Bode plots.

%% System equations in discrete time. Eq. (7.20)
[A,B1,C,D1]=vortex_getsys(Nb,Nw,x-1/4);

% Input matrix with flap hinge at the 3/4 chord.
W(1:Nb,1)=1;
W(1:Nb,2)=2*(xi'-x0);
W(3*Nb/4+1:Nb,3)=1;
W(3*Nb/4+1:Nb,4)=2*(transpose(xi(3*Nb/4+1:Nb))-3/4);

B=B1(:,1:Nb)*W;
D=D1(:,1:Nb)*W;

% Construct state-space model of vortex system, normalize outputs by 2*pi.
sys=ss(A,B,C/2/pi,D/2/pi,2/Nb);


%% Compute discrete-time eigenvalues.
% The discrete eigenvalues are very closed to the unit circle while the 
% continuous eigenvales are very close to the imaginary axis. Yet the 
% system is stable.
[V,lambdafull]=eig(full(A));
for j=1:size(A(:,1))         % Retain only non-zero eigenvalues.
    if abs(lambdafull(j,j)) > 1e-6
    lambda(j)=lambdafull(j,j);
    end
end

% Zero-order hold transformation.
lambdat=log(lambda)/(2/Nb);

% Plot system eigenvaues in continuous time. 
figure(10)
plot(real(lambdat),imag(lambdat),'b.')
xlabel('re(\lambda)')
ylabel('im(\lambda)')
title(['Max root' num2str(max(real(lambdat)))])

%% Compute bode diagram for all i/o in the problem.
[magn,phasen]=bode(sys,wbode);


%% Analytical solutions.
% Before we compute the analytical expressions, define Theodorsen's 
% function as a frequency-response data model object.
k=0:0.0005:10;
thsys=frd(theodorsen(k),k);

iksys=tf([1 0],1);        % Differential operator.


%% Alpha to lift
figure(11)

% Analytical solution (pitch around 1/4 chord).
nu   =-0.5;
theta=pi-acos(nu);

CLqsa0= 2*pi;
CLqsa1= 2*pi*(1/2-nu);
CLnca1= pi;
CLnca2= -nu*pi;

CLa= thsys*(CLqsa0 + CLqsa1*iksys) + (CLnca1*iksys + CLnca2*iksys*iksys);
[magt,phaset]=bode(CLa/2/pi,wbode);

subplot(2,1,1)
 plot(wbode,squeeze(magt),'b-','LineWidth',2), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
subplot(2,1,2)
 plot(wbode,squeeze(phaset),'b-','LineWidth',2), hold on
 xlabel ('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel ('Phase (deg)','FontSize',12,'FontWeight','bold') 

% Solution from vortex method.
Z0=squeeze(magn(1,1,:).*exp(i*phasen(1,1,:)*pi/180));
Z1=squeeze(magn(1,2,:).*exp(i*phasen(1,2,:)*pi/180));
Z= Z0 +i*wbode'.*Z1;

subplot(2,1,1)
 plot(wbode,abs(Z),'k:','LineWidth',2)
subplot(2,1,2)
 plot(wbode,angle(Z)*180/pi,'k:','LineWidth',2)


%% Alpha to moment
figure(12)

% Analytical solution (pitch around 1/4 chord).
nu   =-0.5;
theta=pi-acos(nu);

CMqsa0= 0;
CMqsa1= -(pi/4);
CMnca1= -(pi/4);
CMnca2= -(pi/4)*(1/4-nu);

CMa= CMqsa0 + CMqsa1*iksys + (CMnca1*iksys + CMnca2*iksys*iksys);
[magt,phaset]=bode(CMa/(pi/4),wbode);

subplot(2,1,1)
 plot(wbode,squeeze(magt),'b-','LineWidth',2), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
subplot(2,1,2)
 plot(wbode,squeeze(phaset)-360,'b-','LineWidth',2), hold on
 xlabel ('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel ('Phase (deg)','FontSize',12,'FontWeight','bold') 

% Solution from vortex method.
Z0=squeeze(magn(2,1,:).*exp(i*phasen(2,1,:)*pi/180));
Z1=squeeze(magn(2,2,:).*exp(i*phasen(2,2,:)*pi/180));
Z= Z0 +i*wbode'.*Z1;

subplot(2,1,1)
 plot(wbode,8*abs(Z),'k:','LineWidth',2)
subplot(2,1,2)
 plot(wbode,angle(Z)*180/pi,'k:','LineWidth',2)
 legend('Analytical',['N_b=' int2str(Nb) ', N_w/N_b=' int2str(Nw/Nb)])
 
 
%% Beta to lift.
figure(13)

% Analytical solution (flap hinge at the 3/4 chord).
nu=0.5;
theta=pi-acos(nu);

CLqsb0= 2*pi - 2*theta + 2*sin(theta);
CLqsb1=(2*pi - 2*theta) * (1/2-nu) + sin(theta)*(2-nu);
CLncb1= pi - theta - nu*sin(theta);
CLncb2=-nu*(pi-theta) + (1/3)*(2+nu^2)*sin(theta);

CLb= thsys*(CLqsb0 + CLqsb1*iksys) + (CLncb1*iksys + CLncb2*iksys*iksys);
[magt,phaset]=bode(CLb/2/pi,wbode);

subplot(2,1,1)
 plot(wbode,squeeze(magt),'b-','LineWidth',2), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold') 
subplot(2,1,2)
 plot(wbode,squeeze(phaset),'b-','LineWidth',2), hold on
 xlabel ('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel ('Phase (deg)','FontSize',12,'FontWeight','bold') 

% Solution from vortex method.
Z0=squeeze(magn(1,3,:).*exp(i*phasen(1,3,:)*pi/180));
Z1=squeeze(magn(1,4,:).*exp(i*phasen(1,4,:)*pi/180));
Z= Z0 +i*wbode'.*Z1;

subplot(2,1,1)
 plot(wbode,abs(Z),'k:','LineWidth',2)
subplot(2,1,2)
 plot(wbode,angle(Z)*180/pi,'k:','LineWidth',2)

 
%% Beta to moment.
figure(14)

% Analytical solution (flap hinge at the 3/4 chord).
nu=0.5;
theta=pi-acos(nu);

CMqsb0=-((1+nu)/2)*sin(theta);
CMqsb1=-(1/4)*( pi-theta + (2/3)*(1/2-nu)*(2+nu)*sin(theta));
CMncb1=-(1/4)*( pi-theta + (1/3)*(2-3*nu-2*nu^2)*sin(theta));
CMncb2=-(1/4)*( (1/4-nu)*(pi-theta) ... 
               + (2/3-5*nu/12+nu^2/3+nu^3/6)*sin(theta));

CMb= CMqsb0 + CMqsb1*iksys + (CMncb1*iksys + CMncb2*iksys*iksys);
[magt,phaset]=bode(CMb/(pi/4),wbode);

subplot(2,1,1)
 plot(wbode,squeeze(magt),'b-','LineWidth',2), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
subplot(2,1,2)
 plot(wbode,squeeze(phaset)-360,'b-','LineWidth',2), hold on
 xlabel ('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel ('Phase (deg)','FontSize',12,'FontWeight','bold') 
 
% Solution from vortex method.
Z0=squeeze(magn(2,3,:).*exp(i*phasen(2,3,:)*pi/180));
Z1=squeeze(magn(2,4,:).*exp(i*phasen(2,4,:)*pi/180));
Z= Z0 +i*wbode'.*Z1;

subplot(2,1,1)
 plot(wbode,8*abs(Z),'k:','LineWidth',2)
subplot(2,1,2)
 plot(wbode,angle(Z)*180/pi,'k:','LineWidth',2)
 legend('Analytical',['N_b=' int2str(Nb) ', N_w/N_b=' int2str(Nw/Nb)])
 
% eof