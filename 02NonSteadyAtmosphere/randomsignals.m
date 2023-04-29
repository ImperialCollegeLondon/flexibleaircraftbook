% randomsignals.m 
%  Examples for random processes (using coloured Gaussian noise).
%
%  Notes:
%   - Every time the function is called, it will generate a different 
%     time history.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clear all

%% Correlation function
% We choose first two transfer functions that will be used to shape
% the white noise, and compute their corresponding state-space models.
%
% First TF is von Karman turbulence filter (in the w direction) obtained
% with Campbell's Rational Function Approximation.

G1= tf([sqrt(8/3) 1],[1 1]) ...
  * tf([91/12 52 60],[935/216 561/12 102 60]);
sys1=ss(G1,'minimal');

% Second transfer function is simply 1/(s+10)
G2nodelay= tf([1],[1 10]);
sys2=ss(G2nodelay,'minimal');

% Determine rms for each signal.
sigma_1=0.9;
sigma_2=0.75;

% Obtain time history of Gaussian noise.
Ntimes=10000000;      % Number of time steps (very long, so that the 
                      % the statistics are representative).
dt=0.001;             % Time step. 
t=0:dt:Ntimes*dt;
Nlags=2500;

% Add 0.5 s delay on w2 with respect to w1. Just because we want to see 
% its effect. Remove to see how the results change.
Ndelay=500;
TimeDelay=Ndelay*dt; % Delay on w2

% Simulate time response of the system to the white noise and normalize the
% resulting signal to have unit rms (it is rms(n)=1, but not its PSD. This 
% is why we needed such very large sampling, Ntimes).

n=randn(Ntimes+Ndelay+1,1);  % White noise.
[w1,t1]=lsim(sys1,n(Ndelay+1:Ndelay+Ntimes+1,1),t,'foh');
norm_w1=rms(w1);
w1=sigma_1*w1/norm_w1;

[w2,t2]=lsim(sys2,n(1:Ntimes+1,1),t,'foh'); 
norm_w2=rms(w2);
w2=sigma_2*w2/norm_w2;

% Compute the correlation function (a 2x2 matrix function in this problem)
[r1,lag1]=xcorr(w1,Nlags,'unbiased');
[r2,lag2]=xcorr(w2,Nlags,'unbiased');
[r12,lag12]=xcorr(w1,w2,Nlags,'unbiased');
[r21,lag21]=xcorr(w2,w1,Nlags,'unbiased');

% Plot time histories of both signals only up to t=10.
figure
subplot(2,1,1)
 plot(t1,w1,'k')
 axis([0 10 -4 4])
 ylabel('w_1','FontSize',14,'FontWeight','bold')
subplot(2,1,2)
 plot(t2,w2,'k')
 axis([0 10 -4 4])
 xlabel('t [s]','FontSize',14,'FontWeight','bold')
 ylabel('w_2','FontSize',14,'FontWeight','bold')
 
% Plot the correlation function.
figure 
subplot(2,2,1)
 plot(lag1*dt,r1,'LineWidth',2)
 axis([-2 2 -0.1 1])
 ylabel('\phi_{11}','FontSize',14,'FontWeight','bold')
subplot(2,2,4)
 plot(lag2*dt,r2,'LineWidth',2)
 axis([-2 2 -0.1 1])
 ylabel('\phi_{22}','FontSize',14,'FontWeight','bold')
 xlabel('\tau [s]','FontSize',14,'FontWeight','bold')
subplot(2,2,2)
 plot(lag12*dt,r12,'LineWidth',2)
 axis([-2 2 -0.1 1])
 ylabel('\phi_{12}','FontSize',14,'FontWeight','bold')
subplot(2,2,3)
 plot(lag21*dt,r21,'LineWidth',2)
 axis([-2 2 -0.1 1])
 ylabel('\phi_{21}','FontSize',14,'FontWeight','bold')
 xlabel('\tau [s]','FontSize',14,'FontWeight','bold')


 % Compute the integral time scale for each of the signals.
 T1=trapz(lag1(Nlags+1:end)*dt,r1(Nlags+1:end))
 %subplot(2,2,1), ax = gca; ax.XTick= T1;
 
 T2=trapz(lag2(Nlags+1:end)*dt,r2(Nlags+1:end))
 %subplot(2,2,4), ax = gca; ax.XTick= T2;
 
 % Compute the time delay between both signals (note that this was
 % introduced ``by hand'' with Ndelay).
 [~,I] = max(abs(r12));
 timeDiff=lag12(I)*dt;
 %subplot(2,2,2), ax = gca; ax.XTick= timeDiff;
 %subplot(2,2,3), ax = gca; ax.XTick=-timeDiff;

 
 
%% Spectral Description
% For convenience, we aggregagate both transfer functions to define 
% a system with two inputs and two outputs.
G=tf(parallel(sys1,G2nodelay,1,1,[],[]));

% Compute the PSDs. Since we have the transfer functions, we use
% the spectral factorization theorem. If TFs were not avaialble, we would
% need to carry an fft of the correlation functions above.
PSD(1,1)=G(1,1)*ctranspose(G(1,1));
PSD(1,2)=G(2,1)*ctranspose(G(1,1));
PSD(2,1)=G(1,1)*ctranspose(G(2,1));
PSD(2,2)=G(2,1)*ctranspose(G(2,1));
omega=[0:0.1:20];
[magPSD,phasePSD]=bode(PSD,omega);
PhaseDelay=unwrap(angle(exp(-j*omega/2)))*180/2/pi;  

% Compute norms.
Phi1_norm=trapz(omega,squeeze(magPSD(1,1,:)));
Phi2_norm=trapz(omega,squeeze(magPSD(2,2,:)));


% Plot results.
figure
 subplot(4,2,1)
 plot(omega,sigma_1*sigma_1*squeeze(magPSD(1,1,:))/Phi1_norm,'LineWidth',2)
 title('\Phi_{11}')
 ylabel('mag')
 subplot(4,2,3)
 plot(omega,squeeze(phasePSD(1,1,:)),'LineWidth',2)
 ylabel('phase')
 axis([0 20 -30 30])
 
 subplot(4,2,2)
 plot(omega,sigma_1*sigma_2*squeeze(magPSD(1,2,:))/sqrt(Phi1_norm*Phi2_norm),'LineWidth',2)
 title('\Phi_{12}')
 ylabel('mag')
 subplot(4,2,4)
 plot(omega,squeeze(phasePSD(1,2,:))+PhaseDelay','LineWidth',2)
 ylabel('phase')
 
 subplot(4,2,5)
 plot(omega,sigma_2*sigma_1*squeeze(magPSD(2,1,:))/sqrt(Phi1_norm*Phi2_norm),'LineWidth',2)
 title('\Phi_{21}')
 ylabel('mag')
 subplot(4,2,7)
 plot(omega,squeeze(phasePSD(2,1,:))-PhaseDelay','LineWidth',2)
 ylabel('phase')
 xlabel('\omega [rad/s]')
 
 subplot(4,2,6)
 plot(omega,sigma_2*sigma_2*squeeze(magPSD(2,2,:))/Phi2_norm,'LineWidth',2) 
 ylabel('mag')
 title('\Phi_{22}')

 subplot(4,2,8)
 plot(omega,squeeze(phasePSD(2,2,:)),'LineWidth',2)  
 ylabel('phase')
 axis([0 20 -30 30])
 xlabel('\omega [rad/s]')
 
 % Compute the coherence. It is equal to 1 since both signals came from 
 % the same white noise...
 figure
  for i=1:length(omega)
      coh12(i)=magPSD(1,2,i)^2/(magPSD(1,1,i)*magPSD(2,2,i));
  end
  plot(omega,coh12,'LineWidth',2)
 axis([0 20 0 2])
 % eof