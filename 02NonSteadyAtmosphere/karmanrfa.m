% Function karmanrfa. 
%  Plot analytical expressions for the vertical component of the 
%  velocity in the von Karman turbulence model, and compare with 
%  Campbell's rational function approximation.
%
%  Note: Frequencies are normalized by T=aL/V.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: November 2014. 
%%%%%%%%%%%%%%%%%%%%
close all, clear all

omega=0:0.01:100;  % Nondimensional angular frequency.

% Define frequency-domain analytic function.
fkarman=@(xx) (1 +sqrt(8/3)*xx) ./ (1 + xx).^(11/6);
GKarman  = frd(fkarman(omega),omega);
PSDKarman= (1+(8/3)*(omega).^2) ./ (1+ (omega).^2).^(11/6);

% Campbell's Rational function approximation. Generate PSD values for
% the range of non-dimensional omega defined above.
GCampbell=tf([sqrt(8/3) 1],[1 1]) ...
         *tf([91/12 52 60],[935/216 561/12 102 60]);
PSDCampbell=freqresp(GCampbell*conj(GCampbell),omega);

% Plot PSD. (include abs simply to avoid numerical errors)
figure
loglog(omega,squeeze(abs(PSDKarman)),'b-','LineWidth',2), hold on 
loglog(omega,squeeze(abs(PSDCampbell)),'k--','LineWidth',2)
axis([0.01 100 0.001 2])
legend('von Karman','Campbell approximation')
 ylabel('\Phi_{w3}/C_3^2','FontSize',14,'FontWeight','bold')
 xlabel('\omegaT','FontSize',14,'FontWeight','bold')

% Bode plot of the analytical and Campbell's approximated transfer 
% functions.
figure
[mag,phi]=bode(GKarman,omega);
subplot(2,1,1)
 semilogx(omega,squeeze(mag),'b-','LineWidth',2), hold on
subplot(2,1,2)
 semilogx(omega,squeeze(phi),'b-','LineWidth',2), hold on

[mag,phi]=bode(GCampbell,omega);
subplot(2,1,1)
 semilogx(omega,squeeze(mag),'k--','LineWidth',2)
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 legend('von Karman','Campbell approximation')
subplot(2,1,2)
 semilogx(omega,squeeze(phi),'k--','LineWidth',2)
 axis([0.01 100 -90 5])
 xlabel('\omegaT','FontSize',14,'FontWeight','bold')
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')

% eof