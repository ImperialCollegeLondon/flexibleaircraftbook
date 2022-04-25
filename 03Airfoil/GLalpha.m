% Function GLalpha. 
% Plot transfer function between angle of incidence and lift using 
% unsteady thin aerofoil theory.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: August 2014. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Theodorsen function and function between AoA and lift.
theod=@(xx) besselh(1,2,xx)./(besselh(1,2,xx)+i*besselh(0,2,xx));
GLa=@(xx,nu) 2*pi*theod(xx).*(1+i*xx*(1/2-nu))+pi*(i*xx+nu*xx.^2);

% Plot transfer function
figure 
k=0:0.001:1;
for nu_ea=0:-0.25:-1;
  subplot(2,1,1)
   plot(k,abs(GLa(k,nu_ea)/2/pi),'-','Color',[0 1 1]*abs(nu_ea/1.2), 'LineWidth',2), hold on
  subplot(2,1,2)
   plot(k,phase(GLa(k,nu_ea))*180/pi,'-','Color',[0 1 1]*abs(nu_ea/1.2), 'LineWidth',2), hold on
end

% Write legends.
subplot(2,1,1)
    ylabel('Magnitude (abs)','FontSize',19,'FontWeight','bold')
    grid on
subplot(2,1,2)
    xlabel('Reduced frequency, k','FontSize',19,'FontWeight','bold')
    ylabel('Phase (deg)','FontSize',19,'FontWeight','bold')
    grid on
    
% Plot time-history of alpha vs. CL for varying k.
figure
nu_ea=0.25;
theta=0:0.001:pi/2;  % Full cycle in k*s.
alpha=cos(theta);

k=0.01; plot(alpha,real(GLa(k,nu_ea)*exp(i*theta))/2/pi,'b--','LineWidth',2), hold on
k=0.1;  plot(alpha,real(GLa(k,nu_ea)*exp(i*theta))/2/pi,'b-','LineWidth',2), hold on
k=0.5;  plot(alpha,real(GLa(k,nu_ea)*exp(i*theta))/2/pi,'b:','LineWidth',2)

legend('k=0.01','k=0.1','k=0.5')
xlabel('Angle of incidence, \alpha/\alpha_0','FontSize',16,'FontWeight','bold')
ylabel('Lift coefficient, C_L/(2\pi)','FontSize',16,'FontWeight','bold')
