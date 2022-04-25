% function wagner_kuessner.m
%
% Compute step responses to pitch and gust inputs.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: August 2014. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% Obtain RFA approximation to Theodorsen's function. Remove asymptotic
% value of C(k) as k->infty to speed convergence.
k=0:0.0005:20;
theod=@(x) besselh(1,2,x)./(besselh(1,2,x)+i*besselh(0,2,x))-0.5;
tsys=frd(theod(k),k);
tsys10=fitmagfrd(tsys,10)

figure
[yt10,xt10]=step(0.5+tsys10)


sears=@(x) conj(2./(pi.*x.*(besselh(0,x)+i*besselh(1,x))));
ssys=frd(sears(k),k);
ssys10=fitmagfrd(ssys,10)
figure
[ys10,xs10]= step(ssys10)


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
