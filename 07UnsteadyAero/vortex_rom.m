% vortex_rom.m
%
%  Compute balanced solutions for the state-space vortex description.
%
%  It corresponds to Example 7.2 in Palacios & Cesnik (CUP, 2023)
%   https://doi.org/10.1017/9781108354868
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: August 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

Nb=100;                        % Number of segments on aerofoil
Nw=Nb*30;                      % Number of segments in wake
dx=1/Nb;                       % Nondimensional panel length

% Define non-dimensional coordinates from leading edge.
x=dx/4:dx:1+Nw/Nb-((3*dx)/4);  % Coordinates of vortices (aerofoil/wake).
xi=3*dx/4:dx:1-dx/4;           % Coordinates of the Nb collocation points.

x_ea=0.25;                     % Pitch rotations about 0.25c.
x_fh=0.75;                     % Flap hinge at 0.75c.

wbode=0.001:0.001:2;           % Frequencies for the Bode plots.

%% Aerofoil system equations in discrete time.
% Output is lift and moments about the 1/4 chord (divided by 2*pi).
[Aa,Ba,Ca,Da]=vortex_getsys(Nb,Nw,x-1/4);  
sysa=ss(Aa,Ba,Ca/2/pi,Da/2/pi,2/Nb);

% Contribution of pitch / flap motions to input matrix.
W(1:Nb,1)= 1;
W(1:Nb,2)= 2*(xi'-x_ea)+1/Nb;
for i=1:Nb
    if xi(i) >x_fh
        W(i,3)=1;
        W(i,4)=2*(xi(i)-x_fh)+1/Nb;
    else
        X(i,3:4)=0;
    end
end

% Possibly include gust input to the system.
if 0
    W(1:Nb,5)=0;
    % Define gust system.
    Ag=zeros(Nb,Nb); for j=2:Nb, Ag(j,j-1)=1; end
    Bg=zeros(Nb,5); Bg(1,5)=1;
    Cg=eye(Nb);
    Dg=0;

    sysg=ss(Ag,Bg,Cg,Dg,2/Nb);

    % Combined system.
    sys=sysa*(ss(W)+sysg);
else
    sys=sysa*ss(W);
end


%% Model reduction (balancing)
[sysb, g, T, Ti]=balreal(sys);
% Plot Hankel singular values.
figure
semilogy(g(1:50)/max(g),'ko')
xlabel('Order, j','FontSize',16,'FontWeight','bold')
ylabel('Normalized HSV, \sigma_j','FontSize',16,'FontWeight','bold')

[magsys,phasesys]   =bode(sys,wbode);
[magrsys2,phasesys2]=bode(modred(sysb,3:Nw+Nb),wbode);
[magrsys3,phasesys3]=bode(modred(sysb,4:Nw+Nb),wbode);
[magrsys4,phasesys4]=bode(modred(sysb,5:Nw+Nb),wbode);


%% Plot model reduction individual components.
for k=1:2  % Repeat for each output in the system (lift, moment).
  figure(20+k)
  for j=1:4
    subplot(2,4,j)
     plot(wbode,squeeze(magsys  (k,j,:)),'b-','LineWidth',2), hold on
     plot(wbode,squeeze(magrsys2(k,j,:)),'k-.','LineWidth',2)
     plot(wbode,squeeze(magrsys3(k,j,:)),'k:','LineWidth',2)
     plot(wbode,squeeze(magrsys4(k,j,:)),'k--','LineWidth',2)
     if k==1, axis([0 2 0 1]), end  
     
    subplot(2,4,4+j)
     plot(wbode,squeeze(phasesys (k,j,:)),'b','LineWidth',2), hold on
     plot(wbode,squeeze(phasesys2(k,j,:)),'k-.','LineWidth',2)
     plot(wbode,squeeze(phasesys3(k,j,:)),'k:','LineWidth',2)
     plot(wbode,squeeze(phasesys4(k,j,:)),'k--','LineWidth',2)   
     if k==1 && j<5, axis([0 2 -15 60]), end
     xlabel('k','FontSize',12,'FontWeight','bold')
  end
  % Add legends.
subplot(2,4,1)
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 title('input: \alpha','FontSize',12)
subplot(2,4,2)
 title('input: d\alpha/ds','FontSize',12)
subplot(2,4,3)
 title('input: \delta','FontSize',12)
subplot(2,4,4)
 title('input: d\delta/ds','FontSize',12)
 legend('Full','N=2','N=3','N=4','FontSize',12)
subplot(2,4,5)
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')
end

% eof

