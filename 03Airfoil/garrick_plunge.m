% Function garrick_plunge
%
% Computes induced drag due to harmonic plunging motions.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: May 2018. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all, clear all

% Define frequency-domain analytic functions.
theod=@(xx) besselh(1,2,xx)./(besselh(1,2,xx)+i*besselh(0,2,xx));

linestyles={'b:','b--','b-.','b-'}
k=[0.25:0.25:1];
for i=1:length(k)
    ks=0:0.01:3*pi/2;
    h=cos(ks);
    C=theod(k(i));
    Cd=-2*pi*k(i)^2*(imag(C)*cos(ks)+real(C)*sin(ks)).^2;
    
    plot(h,Cd,linestyles{i},'LineWidth',1), hold on
    grid on
end
ylabel('C_{Di}/(2\pi)')
xlabel('h_0/(c/2)')
legend('k=0.25','k=0.5','k=0.75','k=1')