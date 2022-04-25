% Function karman2d. 
%  Plot normalized von Karman 2-D turbulence models.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2016.
%%%%%%%%%%%%%%%%%%%
close all, clear all

x=0.001:0.001:100;   %omega*T

% Vertical/lateral von Karman turbulence spectra:
Karman_v= (1+(8/3)*x.^2) ./ ((1+ x.^2).^(11/6));
C=2^(1/6)/gamma(5/6);
slist=[0.0001 0.1:0.1:1]';
linecolour=interp1q([0 1]',[0 1]',slist)*0.75;
for ii=1:length(slist)
    s=slist(ii);
    r1=s.*sqrt(1+x.^2);
    r2=s./sqrt(1+x.^2);
    Karman_vv= (8/3*r2.^( 5/6).*besselk( 5/6,r1) ...
               -r2.^(11/6).*besselk(11/6,r1));
               
% Define frequency-domain analytic functions.
  loglog(x,(C*Karman_vv./Karman_v),'Color', linecolour(ii)*[1 1 1], 'LineWidth',2), hold on
  ylabel('Coherence, \gamma_{ww}^{1/2}','FontSize',14,'FontWeight','bold')
  xlabel('\omega T','FontSize',14,'FontWeight','bold')
% legend('von Karman, vertical','Dryden, vertical',...
%        'von Karman, longitudinal', 'Dryden, longitudinal')

axis([0.01 10 0.01 1.1])
end


% Values for omega=0
figure
r=0.0001:0.0001:0.2;
r1=sqrt(2)*r;
r2=r/sqrt(2);
y1=(8/3*r.^(5/6).*besselk(5/6,r)-r.^(11/6).*besselk(11/6,r)).^2;
y2=(8/3*r2.^(5/6).*besselk(5/6,r1)-r2.^(11/6).*besselk(11/6,r1)).^2;
y1=y1/max(y1);
y2=y2/max(y2);
plot(r,y1,'k-.','LineWidth',2), hold on
plot(r,y2,'k--','LineWidth',2)
ylabel('Coherence, \gamma_{ww}','FontSize',14,'FontWeight','bold')
xlabel('\Delta y/al','FontSize',14,'FontWeight','bold')

p1=polyfit(r,y1,3)
plot(r,p1(1)*r.^3+p1(2)*r.^2+p1(3)*r+p1(4),'k-','LineWidth',1)
p2=polyfit(r,y2,3)
plot(r,p2(1)*r.^3+p2(2)*r.^2+p2(3)*r+p2(4),'k-','LineWidth',1)
legend('von Karman, \omega=0','von Karman, \omegaT=1','3^{rd}-order polynomials')

% Davenport approximation

c=8;
figure
for ii=1:length(slist)
    s=slist(ii);
    Davenport_vv=exp(-2*pi*c*s*x);
               
% Define frequency-domain analytic functions.
  loglog(x,Davenport_vv,'Color', linecolour(ii)*[1 1 1], 'LineWidth',2), hold on
  ylabel('Coherence, \gamma_{ww}^{1/2}','FontSize',14,'FontWeight','bold')
  xlabel('\omega T','FontSize',14,'FontWeight','bold')

axis([0.01 10 0.01 1.1])
end
