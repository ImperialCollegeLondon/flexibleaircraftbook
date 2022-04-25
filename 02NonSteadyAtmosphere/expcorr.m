% Function expcorr
% Compute the PSD of a mass/spring to an exponentially-correlated noise.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2022. 
%%%%%%%%%%%%%%%%%%%
clear all, close all

% Damping coefficient.
damp=0.1;

% Relative scale \omega_0\theta
leg={'\omega_0\theta=0.5',...
     '\omega_0\theta=1',...
     '\omega_0\theta=2'};
theta=[0.5 1 2];
lines={'k-','k--','k:'};
 
% FRF response of the system
omega=0:0.001:3;
H=1./(1-omega.^2);
sys=tf([1],[1 2*damp 1]);

[Hamp,Hphase]=bode(sys,omega);

figure
subplot(2,1,1)
plot(omega,squeeze(Hamp),'k','LineWidth',2)
ylabel('magnitude(H)','FontSize',14)
subplot(2,1,2)
plot(omega,squeeze(Hphase),'k','LineWidth',2)
xlabel('\omega/\omega_0','FontSize',14)
ylabel('phase(H) [deg]','FontSize',14)
axis([0 3 -180 0])


% PSD of the input signal

figure
for k=1:length(theta)
    PSDw=1./(1+(theta(k)^2)*omega.^2);

    plot(omega,PSDw,lines{k},'LineWidth',2), hold on
end
xlabel('\omega/\omega_0','FontSize',16)
ylabel('\Phi_w/\Phi_{w,max}','FontSize',16)
legend(leg);


% PSD of the output

figure
for k=1:length(theta)
    PSDy=1./(1+(theta(k)^2)*omega.^2) ...
        .*1./((1-omega.^2).^2+4*damp^2*omega.^2);

    plot(omega,PSDy,lines{k},'LineWidth',2), hold on
end
xlabel('\omega/\omega_0','FontSize',16)
ylabel('\Phi_y/\Phi_{w,max}','FontSize',16)
legend(leg);

% eof