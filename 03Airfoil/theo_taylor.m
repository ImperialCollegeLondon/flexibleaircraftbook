clear all, close all
syms x
f3=series(besselk(1,x)./(besselk(0,x)+besselk(1,x))-1,x,0,'Order',3)
f5=series(besselk(1,x)./(besselk(0,x)+besselk(1,x))-1,x,0,'Order',5)


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
 legend('Theodorsen, C(k)','Taylor 3', 'Taylor 5')

