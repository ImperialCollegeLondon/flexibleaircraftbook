% Function garrick_plunge
%
% Computes induced drag due to harmonic plunging motions.
%
% % Dependencies:
%    theodorsen.m: Analytical expression for Theodorsen's lift deficiency
%                  function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all, clear all

linestyles={'b:','b--','b-.','b-'}
k=[0.25:0.25:1];
for i=1:length(k)
    ks=0:0.01:2*pi;
    h=cos(ks);
    C=theodorsen(k(i));
    Cd=-2*pi*k(i)^2*(imag(C)*cos(ks)+real(C)*sin(ks)).^2;
    
    plot(h,Cd,linestyles{i},'LineWidth',1), hold on
    grid on
end
ylabel('C_{Di}/(2\pi)')
xlabel('h_0/(c/2)')
legend('k=0.25','k=0.5','k=0.75','k=1')