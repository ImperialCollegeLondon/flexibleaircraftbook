% Function karmandryden. 
%  Plot normalized von Karman and Dryden turbulence models.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: November 2014. 
%%%%%%%%%%%%%%%%%%%

close all, clear all

x=0.01:0.0001:100;   %kappa*L (nondimensional wavenumber)

% Vertical/lateral and horizontal von Karman turbulence spectra:
Karman_v= (1+(8/3)*(1.339*x).^2) ./ (1+ (1.339*x).^2).^(11/6);
Karman_h= 2 ./ ((1+ 1.339*x).^2).^(5/6);

% Vertical/lateral and horizontal Dryden turbulence spectra:
Dryden_v= (1+3*x.^2) ./ (1+x.^2).^2;
Dryden_h= 2 ./ (1+x.^2).^2;


% Define frequency-domain analytic functions.
figure
 loglog(x,Karman_v,'b-', 'LineWidth',3), hold on
 loglog(x,Dryden_v,'b-.','LineWidth',3)
 loglog(x,Karman_h,'k-', 'LineWidth',2), hold on
 loglog(x,Dryden_h,'k-.','LineWidth',2)
 ylabel('\Phi^*_w/(\sigma_w^2L/\pi)','FontSize',14,'FontWeight','bold')
 xlabel('\kappa L','FontSize',14,'FontWeight','bold')
 legend('von Karman, vertical','Dryden, vertical',...
        'von Karman, longitudinal', 'Dryden, longitudinal')

axis([0.1 100 0.001 2])