%
% gust_step: Response of an unsupported airfoil to a sharped-edged gust.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all 

% Obtain RFA approximation to Sears
sears=@(x) conj(2./(pi.*x.*(besselh(0,x)+i*besselh(1,x)))).*exp(-i*x);
k=0:0.02:50;
sysse=frd(sears(k),k);
ssys5=fitmagfrd(sysse,5);

% Step gust.
[xs5,ys5]= step(ssys5)

% Jones's approximation.
a_1=0.165;
a_2=0.335;
b_1=0.0455;
b_2=0.3;

LineStyles={'k-','k:','k-.','k--'};
figure, i=0;
for mu=[5 10 50 100]
i=i+1;
A=[-1/(4*mu+1) -2*a_1*b_1/(4*mu+1) -2*a_2*b_2/(4*mu+1);
    -1 -b_1 0;
    -1 0 -b_2];
B=[2/(4*mu+1); 0; 0];
C=[1 0 0];
sysgust=ss(A,B,C,[]);

[y, t]=lsim(sysgust,xs5,ys5);
plot(t(1:end-1),(y(2:end)-y(1:end-1))./(t(2:end)-t(1:end-1)),LineStyles{i},'LineWidth',2)
hold on
end

axis([0 50 0 0.1]), grid on
xlabel('s'), ylabel('d\nu/ds')
legend('\mu=5','\mu=10','\mu=50','\mu=100')