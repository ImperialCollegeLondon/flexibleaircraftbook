clear all, close all
syms x
f3=series(besselk(1,x)./(besselk(0,x)+besselk(1,x))-1,x,0,'Order',2);
f5=series(besselk(1,x)./(besselk(0,x)+besselk(1,x))-1,x,0,'Order',5);


% Bode plots for Theodorsen's functions.
k=0:0.005:0.5;
figure
subplot(2,1,1)
 plot(k,abs(theodorsen(k)),'b-', 'LineWidth',2), hold on
 plot(k,abs(1+subs(f3,x,1i*k)),'b--', 'LineWidth',2), hold on
 plot(k,abs(1+subs(f5,x,1i*k)),'b:', 'LineWidth',2), hold on

 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 axis([0 0.5 0 1])
 grid on

subplot(2,1,2)
 plot(k,angle(theodorsen(k))*180/pi,'b-', 'LineWidth',2), hold on
 plot(k,angle(1+subs(f3,x,1i*k))*180/pi,'b--','LineWidth',2)
 plot(k,angle(1+subs(f5,x,1i*k))*180/pi,'b:','LineWidth',2)
 xlabel('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')
 %axis([0 1 -45 0])
 grid on
 legend('Analytical','3rd-order series', '5th-order series')



% Compute GAFs for 2-DoF airfoil with a flap.
c=1;                  % Not needed, but to make it clearerl
x_ea=0.35*c;          % Elastic axis, from the LE
x_fh=0.85*c;          % Flap hinge location, from the LE
nu_ea=(x_ea-c/2)/(c/2);
nu_fh=(x_fh-c/2)/(c/2);
theta_fh=acos(-nu_fh);

% Aerodynamics influence coefficient matrices 
% (Eq. 3.37 & 3.45)
CLa0qs=2*pi;
CLa1qs=2*pi*(1/2-nu_ea);
CMa1qs=-pi/4;
CLa1nc=pi;
CLa2nc=-pi*nu_ea;
CMa1nc=-pi/4;
CMa2nc=-pi/4*(1/4-nu_ea);

% (Eq. 3.38 & 3.46)
CLd0qs= 2*pi - 2*theta_fh + 2*sin(theta_fh);
CLd1qs= (1/2-nu_fh)*(2*pi-2*theta_fh) ...
      + (2-nu_fh)*sin(theta_fh);
CMd0qs=-(1/2)*(1+nu_fh)*sin(theta_fh);
CMd1qs=-(1/4)*(pi-nu_fh + (2/3)*(1/2-nu_fh)*(2+nu_fh)*sin(theta_fh));
CLd1nc=pi-theta_fh-nu_fh*sin(theta_fh);
CLd2nc=-nu_fh*(pi-theta_fh) + (1/3)*(2+nu_fh^2)*sin(theta_fh);
CMd1nc=-(1/4)*(pi-theta_fh+(2/3-nu_fh-(2/3)*nu_fh^2)*sin(theta_fh));
CMd2nc=-(1/4)*((1/4-nu_fh)*(pi-theta_fh) ...
              +(2/3-(5/12)*nu_fh+(1/3)*nu_fh^2+(1/6)*nu_fh^3) ...
               *sin(theta_fh));

% Eq. (3.60)          
A_0 = [[0 CLa0qs CLd0qs];
       [0 0      CMd0qs]];
A_1 = [[CLa0qs CLa1qs+CLa1nc CLd1qs+CLd1nc];
       [0      CMa1qs+CMa1nc CMd1qs+CMd1nc]];
A_2 = [[CLa1nc CLa2nc CLd2nc];
       [CMa1nc CMa2nc CMd2nc]];
A_3 = [[0 CLa1qs CLd1qs];
       [0      0      0]];
A_4 = [[CLa0qs CLa1qs CLd1qs];
       [0      0      0]];

% Build a sample of matrices
ksmall=1e-6;
klarge=5;
k=[0 ksmall 0.02:0.02:4.98 klarge-ksmall klarge];
for j=1:length(k)
    ik=1i*k(j);
    A(j,:,:)=A_0+ik*A_1+ik^2*A_2+(theodorsen(k(j))-1)*(A_3+ik*A_4);
end

figure
subplot(2,3,1)
  yyaxis left,  plot(k,squeeze(abs(A(:,1,1))),'b'), hold on
  title('A_{11}')
  yyaxis right, plot(k,squeeze(angle(A(:,1,1)))*180/pi,'r')
   ylim([0 180]), yticks([0 90 180]) 
subplot(2,3,2)
  yyaxis left,  plot(k,squeeze(abs(A(:,1,2))),'b'), hold on
  title('A_{12}')
  yyaxis right, plot(k,squeeze(angle(A(:,1,2)))*180/pi,'r')
  ylim([0 180]), yticks([0 90 180]) 
subplot(2,3,3)
  yyaxis left,  plot(k,squeeze(abs(A(:,1,3))),'b'), hold on
  title('A_{13}')
  yyaxis right, plot(k,squeeze(angle(A(:,1,3)))*180/pi,'r')
  ylim([0 180]), yticks([0 90 180]), legend('Magnitude','Phase (deg)')
subplot(2,3,4)
  yyaxis left,  plot(k,squeeze(abs(A(:,2,1))),'b'), hold on
  title('A_{21}')
  yyaxis right, plot(k,squeeze(angle(A(:,2,1)))*180/pi,'r')
   ylim([-180 0]), yticks([-180 -90 0]) 
subplot(2,3,5)
  yyaxis left,  plot(k,squeeze(abs(A(:,2,2))),'b'), hold on
  title('A_{22}')
  yyaxis right, plot(k,squeeze(angle(A(:,2,2)))*180/pi,'r')
   ylim([-180 0]), yticks([-180 -90 0]) 
subplot(2,3,6)
  yyaxis left,  plot(k,squeeze(abs(A(:,2,3))),'b'), hold on
  title('A_{22}')
  yyaxis right, plot(k,squeeze(angle(A(:,2,3)))*180/pi,'r')
   ylim([-180 0]), yticks([-180 -90 0]) 


% Identified functions
A_0i=reshape(A(1,:,:),2,3);
A_1i=real(squeeze(A(end,:,:)-A(end-1,:,:))/(1i*ksmall))
%temp=(ik*A_1+ik^2*A_2+(theodorsen(ksmall)-1)*(A_3+ik*A_4))/ik;
%A_1i=temp-double(eulergamma)*log(ik/2)*A_3
klarge=10; ik=1i*klarge;
A_2i=real(squeeze(A(end,:,:))-A_0i-1i*k(end)*squeeze(A_1i(end,:,:))/(-k(end)^2));
% eof   