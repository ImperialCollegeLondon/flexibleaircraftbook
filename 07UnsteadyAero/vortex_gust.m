% Function vortex_gust. 
% Compare transfer functions between leading edge gust velocity and 
% lift as obtained by the state-space vortex description and the 
% analytical formulas by Sears.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: July 2015. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
Nb=20;                         % Number of panels on aerofoil
Nw=Nb*30;                      % Number of chordlengths in wake
dx=1/Nb;                       % Nondimensional panel length.
x=dx/4:dx:1+Nw/Nb-((3*dx)/4);  % Coordinates of vortices (aerofoil/wake).
xi=3*dx/4:dx:1-dx/4;           % Coordinates of collocation points.

wbode=0.001:0.001:2;           % Frequencies for the Bode plots.

% System equations for the aerofoil system.
[A,B,C,D]=vortex_getsys(Nb,Nw,x-1/4);
sysa=ss(A,B,C/2/pi,D/2/pi,2/Nb);

% Define gust system.
Ag=zeros(Nb,Nb); for j=2:Nb, Ag(j,j-1)=1; end
Bg=zeros(Nb,1); Bg(1)=1;
Cg=eye(Nb);
Dg=0;

sysg=ss(Ag,Bg,Cg,Dg,2/Nb);

% Combined system.
sysag=sysa*sysg;

%% Bode plots
% Sears function.
k=0:0.0005:20;
sears=@(x) conj(2./(pi.*x.*(besselh(0,x)+i*besselh(1,x)))).*exp(-i*x);
sysse=frd(sears(k),k);
%sysse10=fitmagfrd(sysse,10)

[mag,pag]=bode(sysag,wbode);
[mse,pse]=bode(sysse,wbode);

figure
 subplot(2,1,1)
 plot(wbode,squeeze(mse),'b-','LineWidth',2), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 axis([0 2 0 1])
 
 subplot(2,1,2)
 plot(wbode,squeeze(pse),'b-','LineWidth',2), hold on
 xlabel ('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel ('Phase (deg)','FontSize',12,'FontWeight','bold') 
 %axis([0 2 0 90])

 subplot(2,1,1)
 plot(wbode,squeeze(mag(1,1,:)),'k:','LineWidth',2)
 subplot(2,1,2)
 plot(wbode,squeeze(pag(1,1,:)),'k:','LineWidth',2)
 
 legend('Analytical','N_b=20')

% eof
 
%
