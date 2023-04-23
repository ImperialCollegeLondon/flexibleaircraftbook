% Function theosears. 
%  Plot analytical expressions for Theodorsen's and Sear's functions, 
%  and compute and compare rational function approximations to both 
%  functions using fitmagfrd.
%
% Dependencies:
%    theodorsen.m: Analytical expression for Theodorsen's lift deficiency
%                  function.
%    sears.m: Analytical expression for Sears's function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all, clear all

% Redefine Theodorsen's function
theod=@(xx) theodorsen(xx);

% Bode plots.
k=0:0.0005:3;
figure
subplot(2,1,1)
 plot(k,abs(theod(k)),'b-', 'LineWidth',2), hold on
 plot(k,abs(sears(k)),'k--','LineWidth',2)
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 axis([0 3 0 1])
 grid on
 
subplot(2,1,2)
 plot(k,angle(theod(k))*180/pi,'b-', 'LineWidth',2), hold on
 plot(k,angle(sears(k))*180/pi,'k--','LineWidth',2)
 xlabel('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')
 axis([0 3 -45 0])
 grid on
 legend('Theodorsen, C(k)','Sears, S_0(k)')


%% Rational function approximations to Theodorsen.
deltak=0.01;
k=0:deltak:5;
Wk=ones(size(k));
Wk(1:1/deltak)=100;
systh=frd(theod(k)-0.5,k);
systh2=fitmagfrd(systh,2,1,Wk)+0.5;
systh4=fitmagfrd(systh,4,1,Wk)+0.5;
jones=1-tf([0.165 0], [1 0.0455])-tf([0.335 0], [1 0.3]);

% Bode plots.
k=0:0.005:2;
figure
subplot(2,1,1)
 plot(k,abs(theod(k)),'b-', 'LineWidth',1), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 axis([0 2 0 1])
subplot(2,1,2)
 plot(k,angle(theod(k))*180/pi,'b-', 'LineWidth',1), hold on
 xlabel('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')
 axis([0 2 -20 0])
 
[mag,phi]=bode(jones,k);
subplot(2,1,1)
 plot(k,squeeze(mag),'k:','LineWidth',1)
subplot(2,1,2)
 plot(k,squeeze(phi),'k:','LineWidth',1)

[mag,phi]=bode(systh2,k);
subplot(2,1,1)
 plot(k,squeeze(mag),'k-.','LineWidth',1)
subplot(2,1,2)
 plot(k,squeeze(phi),'k-.','LineWidth',1)
 
k=0:0.05:2;
[mag,phi]=bode(systh4,k);
subplot(2,1,1)
 plot(k,squeeze(mag),'ks','LineWidth',1,'MarkerSize',3)
subplot(2,1,2)
 plot(k,squeeze(phi),'ks','LineWidth',1,'MarkerSize',3)
 
legend('Theodorsen', 'Jones RFA', '2-state RFA', '4-state RFA')
grid on

%% Rational function approximations to Sears.
deltak=0.02;
k=0:deltak:40;  % Note that a much larger value of k_max is needed to 
                % fit the asymptotic value through fitmagfrd (which does
                % enforce that S goes to zero automatically).Wk=ones(size(k));
Wk=ones(size(k));
Wk(1:1/deltak)=100;

sysse=frd(sears(k),k);
sysse2=fitmagfrd(sysse,2,[],Wk);
sysse4=fitmagfrd(sysse,4,[]);
jones=1-tf([0.5 0], [1 0.13])-tf([0.5 0], [1 1]);

% Bode plots.
k=0:0.005:2;
figure
subplot(2,1,1)
 plot(k,abs(sears(k)),'b-', 'LineWidth',1), hold on
 ylabel('Magnitude (abs)','FontSize',12,'FontWeight','bold')
 axis([0 2 0 1]), grid on
subplot(2,1,2)
 plot(k,angle(sears(k))*180/pi,'b-', 'LineWidth',1), hold on
 xlabel('Reduced frequency, k','FontSize',12,'FontWeight','bold')
 ylabel('Phase (deg)','FontSize',12,'FontWeight','bold')
 axis([0 2 -45 0]), grid on
 
[mag,phi]=bode(jones,k);
subplot(2,1,1)
 plot(k,squeeze(mag),'k:','LineWidth',1)
subplot(2,1,2)
 plot(k,squeeze(phi),'k:','LineWidth',1)

[mag,phi]=bode(sysse2,k);
subplot(2,1,1)
 plot(k,squeeze(mag),'k-.','LineWidth',1)
subplot(2,1,2)
 plot(k,squeeze(phi),'k-.','LineWidth',1)

k=0:0.05:2;
[mag,phi]=bode(sysse4,k);
subplot(2,1,1)
 plot(k,squeeze(mag),'ks','LineWidth',1,'MarkerSize',3)
subplot(2,1,2)
 plot(k,squeeze(phi),'ks','LineWidth',1,'MarkerSize',3)
 
legend('Sears', 'Jones RFA', '2-state RFA', '4-state RFA')
