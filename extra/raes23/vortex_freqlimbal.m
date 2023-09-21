% vortex_freqlimbal.m
%
%  Compute FL balanced solutions for the state-space vortex description.
%
%  It has been modified from Example 7.2 in Palacios & Cesnik (CUP, 2023)
%   https://doi.org/10.1017/9781108354868
%
%  A single dependency on each of the 3 degree of freedom is now included, 
%  with derivatives introduce through a second-oder operator.
%
%  This results appeared in a RAeS'23 conference paper
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: Sept 2023. Needs Matlab2023b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
addpath('../../07UnsteadyAero')

Nb=50;                         % Number of segments on aerofoil
Nw=Nb*30;                      % Number of segments in wake 
dx=1/Nb;                       % Nondimensional panel length
ds=2*dx;                       % Nondimensional time step.
knyq=pi/ds;                    % Nyquist reduced frequency.


% Define non-dimensional coordinates from leading edge.
x=dx/4:dx:1+Nw/Nb-((3*dx)/4);  % Coordinates of vortices (aerofoil/wake).
xi=3*dx/4:dx:1-dx/4;           % Coordinates of the Nb collocation points.
                    
x_ea=0.35;                     % Pitch rotations about 0.35c.
x_fh=0.85;                     % Flap hinge at 0.85c.
 
wbode=0.05:0.05:10;            % Frequencies for the Bode plots.

Nrom=[4 6 8];
linestyles={'--',':','-.'};

%% Aerofoil system equations in discrete time.
[Aa,Ba,Ca,Da]=vortex_getsys(Nb,Nw,x-1/4);
sysa=ss(Aa,Ba,Ca/2/pi,Da/2/pi,ds);

% Contribution of airfoil motions to input matrix.
z=tf('z',ds);
D1=1/(2*ds)*(1/z^2)*(3*z^2-4*z+1);

W(1:Nb,1)= ones(Nb,1)+D1*(2*(xi'-x_ea));
W(1:Nb,2)= D1*ones(Nb,1);
for j=1:Nb
    if xi(j) >x_fh
        W(j,3)=1+D1*(2*(xi(j)-x_fh));
    else
        W(j,3)=0;
    end
end

sys=sysa*W;

% Obtain FRF for reference system.
[magsys,phasesys]   =bode(sys,wbode);


figure(11)
for kk=1:2  % Repeat for each output in the system (lift, moment).
  for jj=1:3 
    subplot(2,3,(kk-1)*3+jj)
       Z=squeeze(magsys(kk,jj,:).*exp(1i*phasesys(kk,jj,:)*pi/180));
       yyaxis left,  plot(wbode,real(Z),'-','LineWidth',2)
       xlabel('k','FontSize',12,'FontWeight','bold'), hold on
       ylabel(['A' int2str(kk) int2str(jj)],'FontSize',12,'FontWeight','bold')
       yyaxis right, plot(wbode,imag(Z),'-'), hold on
  end
end


%% Frequency-limited balancing
R2 = reducespec(sys,"balanced");
R2.Options.FreqIntervals = [0,5];
R2 = process(R2);

% Plot FRF for 3 different ROM sizes.
figure(11)
for ii=1:length(Nrom)
    rsys = getrom(R2,Order=Nrom(ii),Method="matchDC");
    [magrsys,phasersys]   =bode(rsys,wbode);
    for kk=1:2  % Repeat for each output in the system (lift, moment).
      for jj=1:3 
        subplot(2,3,(kk-1)*3+jj)
           Z=squeeze(magrsys(kk,jj,:).*exp(1i*phasersys(kk,jj,:)*pi/180));
           yyaxis left,  plot(wbode,real(Z),linestyles{ii},'LineWidth',2)
           xlabel('k','FontSize',12,'FontWeight','bold')
           ylabel(['A' int2str(kk) int2str(jj)],'FontSize',12,'FontWeight','bold')
           yyaxis right, plot(wbode,imag(Z),linestyles{ii}','LineWidth',2)
      end
    end
end
legend('Full',[int2str(Nrom(1)) ' states'],...
              [int2str(Nrom(2)) ' states'],...
              [int2str(Nrom(3)) ' states'])


% eof

