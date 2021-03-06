% Function vortex_rom. 
% Compute balanced solutions for the state-space vortex description.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: July 2015. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

Nb=100;                        % Number of panels on aerofoil
Nw=Nb*30;                      % Number of chordlengths in wake
dx=1/Nb;                       % Nondimensional panel length.
x=dx/4:dx:1+Nw/Nb-((3*dx)/4);  % Coordinates of vortices (aerofoil/wake).
xi=3*dx/4:dx:1-dx/4;           % Coordinates of the Nb collocation points.

x0=1/4;                        % Rotate about the quarter chord.

wbode=0.001:0.001:2;           % Frequencies for the Bode plots.

%% Aerofoil system equations.
[Aa,Ba,Ca,Da]=vortex_getsys(Nb,Nw,x-1/4);
sysa=ss(Aa,Ba,Ca/2/pi,Da/2/pi,2/Nb);

% Contribution of airfoil motions to input matrix.
W(1:Nb,1)=1;
W(1:Nb,2)=2*(xi'+x0-1/2)+1/Nb; % Check whether +1/(Nb) needs be removed.
W(3*Nb/4+1:Nb,3)=1;
W(3*Nb/4+1:Nb,4)=2*(transpose(xi(3*Nb/4+1:Nb))+1/(2*Nb)-3/4);

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
%[magrsys8,phaserys8]=bode(modred(sysb,9:Nw+Nb),wbode);

% Fix results manually...
%phasesys3(1,5,:)=phasesys3(1,5,:)-360;
%phasesys4(1,5,:)=phasesys4(1,5,:)-360;



%% Plot model reduction individual components.
for k=1:2  % Repeat operation for each output in the system (lift, moment).
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
 title('input: \beta','FontSize',12)
subplot(2,4,4)
 title('input: d\beta/ds','FontSize',12)
 legend('Full','N=2','N=3','N=4','FontSize',12)
subplot(2,4,5)
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')
end

% eof

