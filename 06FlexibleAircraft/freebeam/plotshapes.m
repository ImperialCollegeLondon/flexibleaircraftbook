%
% plotshapes.m: Plot mode shapes.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all, close all
L=1;
xrange=0:L/200:L;
linestyles={'k-','k--','k-.','k:'};
for j=1:4
  Nfun = @(x)shape(x,L,j);
  plot(xrange,Nfun(xrange),linestyles{j},'LineWidth',2), hold on
  xlabel('\Deltax/(\Delta L)','FontSize',14)
  ylabel('N','FontSize',14)
  grid on
end
legend('N_1','N_2','N_3','N_4')