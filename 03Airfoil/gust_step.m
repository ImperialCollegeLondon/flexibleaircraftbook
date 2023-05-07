% gust_step.m
%  
%    Response of an unsupported airfoil to a sharped-edged gust.
%
% Dependencies:
%    theodorsen_ajbj.m: aj & bj coefficients of RFA to Theodorsen.
%    sears_rfa.m: RFA to Sears's function.
%
% Copyright, Rafael Palacios, May 2023
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all 

% Obtain RFA approximation to Sears
ssys5=sears_rfa(5);

% Step gusts (Kuessner function)
Tmax=100;
[xs5,ys5]= step(ssys5,Tmax)

% RFA of Theodorsen's Lift deficiency function.
[a,b]=theodorsen_ajbj(2);
a_1=a(1); a_2=a(2); b_1=b(1);b_2=b(2);

% Uncomment to use Jones's coefficients instead.
% a_1=0.165; a_2=0.335; b_1=0.0455; b_2=0.3;


% Select parametrization.
LineStyles={'k-','k:','k-.','k--'};
mulist=[5 20 50 100];
i=0;

% Run through all parameters.
for mu=mulist
    i=i+1;
    % Define state-space matrices for current value of mu.
    A=[-1/(4*mu+1) 2*a_1*b_1/(4*mu+1)  2*a_2*b_2/(4*mu+1);
        -1 -b_1 0;
        -1 0 -b_2];
    B=[2/(4*mu+1); 0; 0];
    C=[1 0 0];
    sysgust=ss(A,B,C,[]);

    % Compute the time-history of the velocities
    [y, t]=lsim(sysgust,xs5,ys5);
    figure(1)
    plot(t,y,LineStyles{i},'LineWidth',2)
    hold on 

    % Compute the time-history of the normalized accelerations
    ydot=diff(y)*length(t)/Tmax;
    figure(2)
    plot(t(1:end-1),2*mu*ydot,LineStyles{i},'LineWidth',2)
    hold on 
end

% Add labels to the figures.
figure (1), ylabel('nondim velocity, \nu')
figure (2), ylabel('gust alleviation factor, 2\mu d\nu/ds')

for i=1:2
    figure(i),  xlabel('reduced frequency, s'),
    axis([0 Tmax 0 1]), grid on
    legend(['\mu=' num2str(mulist(1))], ...
           ['\mu=' num2str(mulist(2))], ...
           ['\mu=' num2str(mulist(3))], ...
           ['\mu=' num2str(mulist(4))])
end
