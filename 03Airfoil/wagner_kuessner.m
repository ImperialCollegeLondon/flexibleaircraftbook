% function wagner_kuessner.m
%
% Compute step responses to pitch and gust inputs.
%
% Dependencies:
%    theodorsen.m: Analytical expression for Theodorsen's lift deficiency
%                  function.
%    sears.m: Analytical expression for Sears's function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Obtain RFA approximation to Theodorsen's function. Remove asymptotic
% value of C(k) as k->infty to speed up convergence.
k=0:0.0005:20;
theodm=@(x) theodorsen(x)-0.5;
tsys=frd(theodm(k),k);
tsys10=fitmagfrd(tsys,10)

[yt10,xt10]=step(0.5+tsys10);

ssys=frd(sears(k),k);
ssys10=fitmagfrd(ssys,10)
[ys10,xs10]= step(ssys10);


% Plot results, including Jones approximations
figure
plot(xt10,yt10,'k-', 'LineWidth',2), hold on
plot(xt10,1-0.165*exp(-0.0455*xt10)-0.335*exp(-0.3*xt10),'b-.')
plot(xs10,ys10,'k--','LineWidth',2)
plot(xt10,1-0.5*exp(-0.13*xt10)-0.5*exp(-1*xt10),'b:')
axis([0 50 0 1]); grid on
legend('Wagner, \phi','Wagner, 2-state','Kuessner, \psi','Kuessner, 2-state')
xlabel('s=2V_\infty t/c')
ylabel('\phi(s), \psi(s)')

% eof
