% Function static_flap 
%  Compute ideal lift and moment coefficients from trailing-edge flap.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2022.
%%%%%%%%%%%%%%%%%%%
clear all, close all

% Nondimensional flap-hinge location:
theta_fh=0:0.001:pi;
nu_fh=-cos(theta_fh);

CL_beta0= 2*(pi-theta_fh+sin(theta_fh));
CL_beta1= 2*(pi-theta_fh).*(1/2-nu_fh)+(2-nu_fh).*sin(theta_fh);
CM_beta0=-((1+nu_fh)/2).*sin(theta_fh);
CM_beta1=-(1/4)*(pi-theta_fh+(2/3-nu_fh-2/3*nu_fh.^2).*sin(theta_fh));


[haxes,hline1,hline2] = plotyy(nu_fh,CL_beta1/pi,[nu_fh',nu_fh'],[(4/pi)*CM_beta1',(4/pi)*CM_beta0'])
hold on
plot(nu_fh,CL_beta0/pi,'LineStyle','-','LineWidth',2)
set(hline2,'LineWidth',2);
set(hline1,'LineStyle','--','LineWidth',2);
ylabel(haxes(1),'C_L/\pi'     ,'FontSize',12,'FontWeight','bold')      % label left y-axis
ylabel(haxes(2),'C_M/(\pi/4)' ,'FontSize',12,'FontWeight','bold')      % label right y-axis
xlabel(haxes(2),'\nu_{fh}','FontSize',12,'FontWeight','bold') % label x-axis
legend('C_{L\beta`}','C_{L\beta}','C_{M\beta`}','C_{M\beta}')

set(haxes(1),'YTick',[0 1 2 3 4])
set(haxes(2),'YTick',[-1 -0.75 -0.5 -0.25 0])
%grid on