%
% pitchpsd.m: Plot the coefficients for the pitch response problem.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
styles={'b-','b--','b-.','b:'};
kt=[1 5 10 100];
for k=1:length(kt)
    xx=0:0.01:10;
    phi=1./(1+(xx/kt(k)).^2).* ...
        (1+sqrt(8/3)*xx.^2) ./ (1 + xx.^2).^(11/6);
    plot(xx,phi,styles{k},'LineWidth',2'), hold on
end

legend(['k_t=' num2str(kt(1))],...
       ['k_t=' num2str(kt(2))],...
       ['k_t=' num2str(kt(3))],...
       ['k_t=' num2str(kt(4))])
   
 xlabel('\omega T','FontSize',14,'FontWeight','bold')
 ylabel('\Xi','FontSize',14,'FontWeight','bold')