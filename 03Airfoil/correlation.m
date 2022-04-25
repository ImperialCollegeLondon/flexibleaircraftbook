%
% correlation: Plot bivariate correlation function.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all, close all
y1=linspace(-1,1);
y2=linspace(-1,1);
[Y1,Y2] = meshgrid(y1,y2);
k=0; for r=[0 0.5 0.99]; k=k+1; 
 Z = (1/(2*pi*sqrt(1-r^2))) ...
   * exp((-1/(2*(1-r^2)))*(Y1.^2+Y2.^2-2*r*Y1.*Y2));
 subplot(1,3,k)
 contour(Y1,Y2,Z,'k')
 xlabel('y_1/\sigma_1','FontSize',12,'FontWeight','Bold')
 title(['\rho=' num2str(r)],'FontSize',12,'FontWeight','Bold')
 axis equal
end
subplot(1,3,1)
ylabel('y_2/\sigma_2','FontSize',12,'FontWeight','Bold')
% eof