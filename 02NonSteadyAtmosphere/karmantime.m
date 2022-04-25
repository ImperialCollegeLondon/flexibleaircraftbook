% Function karmantime. 
%  Create a time history for the vertical velocity component of the 
%  von Karman turbulence model, from a rational-function approximation.
%
%  Notes:
%   - Frequencies are normalized with the period T =aL/V.
%   - Every time the function is called, it generates a different time 
%      history.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: November 2014. 
%%%%%%%%%%%%%%%%%%%

close all, clear all

V=100;                % Airspeed, m/s.

a=1.339;              % Non-dimensional parameter in von Karman model.
L=2500*0.3048;        % Turbulence scale, 2500 ft.
T=a*L/V;              % Characteristic time scale of interactions.

% Campbell Rational function approximation, in terms of omega*T.
GCampbell=sqrt(T/(pi*a)) ...
         *tf([sqrt(8/3) 1],[1 1]) ...
         *tf([91/12 52 60],[935/216 561/12 102 60]);

% Create state-space (A,B,C,D) system description, in non-dimensional 
% time: t/T.
sysCampbell=ss(GCampbell,'minimal');


% Obtain time history of Gaussian noise.
Ntimes=1000000;       % Number of time steps (very long, so that the 
                      % the statistics are representative).
n=randn(Ntimes+1,1);  % White noise.
dt=0.001;             % Nondimensional time step (dt/T). 
t=0:dt:Ntimes*dt;

% Simulate time response of the system to the white noise and normalize the
% resulting signal to have unit rms (it is rms(n)=1, but not its PSD. This 
% is why we needed such very large sampling, Ntimes). 
[w,t]=lsim(sysCampbell,n,t,'foh'); 
w=w/rms(w);

% Plot only up to nondimensional time 10 of the results.
figure
subplot(2,1,1)
 plot(t,n)
 axis([0 10 -4 4])
 ylabel('n(t)','FontSize',12,'FontWeight','bold')
 subplot(2,1,2)
 plot(t,w,'k')
 axis([0 10 -4 4])
 xlabel('t/T','FontSize',14,'FontWeight','bold')
 ylabel('w(t)/\sigma_w','FontSize',12,'FontWeight','bold')

 % eof