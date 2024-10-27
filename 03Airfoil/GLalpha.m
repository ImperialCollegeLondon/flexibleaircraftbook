% Function GLalpha. 
% Plot transfer function between angle of incidence and lift using 
% unsteady thin aerofoil theory.
%
% Dependencies:
%    theodorsen.m: Analytical expression for Theodorsen's lift deficiency
%                  function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023.
%
% Note: A live script version is also available, see GLalpha_live.mlx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Define transfer function between AoA and lift.
GLa_qs=@(xx,nu) 2*pi*(1+1i*xx*(1/2-nu));
GLa_c=@(xx,nu) theodorsen(xx).*GLa_qs(xx,nu);
GLa_nc=@(xx,nu) pi*(1i*xx+nu*xx.^2);
GLa =@(xx,nu) GLa_c(xx,nu)+GLa_nc(xx,nu);

% Plot transfer function for varying nu_ea.
figure 
k=0:0.001:1;
for nu_ea=0:-0.25:-1;
  subplot(2,1,1)
   plot(k,abs(GLa(k,nu_ea)/2/pi),'-','Color',[0 1 1]*abs(nu_ea/1.2), 'LineWidth',2), hold on
  subplot(2,1,2)
   plot(k,angle(GLa(k,nu_ea))*180/pi,'-','Color',[0 1 1]*abs(nu_ea/1.2), 'LineWidth',2), hold on
end

% Write legends.
subplot(2,1,1)
   ylabel('Magnitude (abs)','Interpreter','latex','FontSize',19)
subplot(2,1,2)
   xlabel('Reduced frequency, $k$','Interpreter','latex','FontSize',19)
   ylabel('Phase (deg)','Interpreter','latex','FontSize',19)

% Plot contributions to the transfer function for nu_ea=-0.25.
figure
k=0.001:0.001:1;
nu_ea=-0.25;

subplot(2,1,1)
   plot(k,abs(GLa_qs(k,nu_ea)/2/pi),'k-.', 'LineWidth',2), hold on
   plot(k,abs(GLa_c(k,nu_ea)/2/pi),'k--', 'LineWidth',2)
   plot(k,abs(GLa_nc(k,nu_ea)/2/pi),'k:', 'LineWidth',2)
   plot(k,abs(GLa(k,nu_ea)/2/pi),'k-', 'LineWidth',3)
   ylabel('Magnitude (abs)','Interpreter','latex','FontSize',19)

subplot(2,1,2)
   plot(k,angle(GLa_qs(k,nu_ea))*180/pi,'k-.', 'LineWidth',2), hold on
   plot(k,angle(GLa_c(k,nu_ea))*180/pi,'k--', 'LineWidth',2)
   plot(k,angle(GLa_nc(k,nu_ea))*180/pi,'k:', 'LineWidth',2)
   plot(k,angle(GLa(k,nu_ea))*180/pi,'k-', 'LineWidth',3)
   xlabel('Reduced frequency, $k$','Interpreter','latex','FontSize',19)
   ylabel('Phase (deg)','Interpreter','latex','FontSize',19)
   yticks([0 30 60 90])
   legend('quasi-steady','circulatory','non-circulatory','total')


% Plot time-history of alpha vs. CL for varying k.
figure
nu_ea=0.25;
theta=0:0.001:2*pi;  % Full cycle in k*s.
alpha=cos(theta);

k=0.01; plot(alpha,real(GLa(k,nu_ea)*exp(i*theta))/2/pi,'k--','LineWidth',2), hold on
k=0.1;  plot(alpha,real(GLa(k,nu_ea)*exp(i*theta))/2/pi,'k-','LineWidth',2), hold on
k=0.5;  plot(alpha,real(GLa(k,nu_ea)*exp(i*theta))/2/pi,'k:','LineWidth',2)

legend('k=0.01','k=0.1','k=0.5')
xlabel('Angle of incidence, $\frac{\alpha}{\alpha_0}$','Interpreter','latex','FontSize',19)
ylabel('Lift coefficient, $\frac{C_L}{2\pi}$','Interpreter','latex','FontSize',19)
