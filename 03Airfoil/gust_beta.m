%
% gust_beta: Compute gust alleviation factor for 1-DoF airfoil.
%
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all 

% Jones's approximation.
a_1=0.165;
a_2=0.335;
b_1=0.0455;
b_2=0.3;

% Obtain RFA approximation to Sears
sears=@(x) conj(2./(pi.*x.*(besselh(0,x)+i*besselh(1,x)))).*exp(-i*x);
k=0:0.02:50;
sysse=frd(sears(k),k);
ssys5=fitmagfrd(sysse,5);

Hlist=[0.4:0.4:100];
mulist=[.2:.2:100];

for j=1:length(Hlist);  % Gradient in semichords.
  H=Hlist(j);
  for k=1:length(mulist);
  mu=mulist(k);
  
  % Create 1-cos gust velocity profile.
  t=0:H/100:6*H;
  wg(1:200)=0.5*(1-cos(pi*t(1:200)/H));
  wg(201:601)=0;
 
  % Obtain lift due to gust velocity.
  [CLg, tg]=lsim(ssys5,wg,t);  % Lift divided by 2\pi.

  % Equations of motion for the 1-dof aerofoil. 
  A_1=[[4*mu+1 2*a_1 2*a_2];[0 1 0];[0 0 1]];
  A_2=[-2 0 0; -1 -b_1 0; -1 0 -b_2];
  A=inv(A_1)*A_2;
  B0=[1; 0; 0];
  B=inv(A_1)*B0;
  
  % Output is the acceleration.
  C=A(1,:);
  D=B(1,1);
  sysgust=ss(A,B,C,D);

  % Solve problem.
  [ydot, t]=lsim(sysgust,CLg,tg);
  beta(j,k)=4*mu*max(ydot);
  end
end

figure
plot(Hlist,beta(:,25:25:500),'k'), hold on
xlabel('Nondimensional gust gradient, 2H/c','FontSize',14,'FontWeight','bold')
ylabel('Gust alleviation factor, \beta','FontSize',14,'FontWeight','bold')
axis([0 50 0 1])

Href=0.01:0.01:100;
for k=1:length(mulist)
 [Kg(k) iKg(k)]=max(spline(Hlist,beta(:,k),Href));
end
hold on
plot(iKg*0.01,Kg,'b--','LineWidth',2)


figure
Kgref=0.88*mulist./(5.3+mulist);
plot(mulist,Kg,'k','LineWidth',2), hold on
plot(mulist,Kgref,'b--','LineWidth',2)
grid on
axis([0 100 0 1])
xlabel('Mass parameter, \mu','FontSize',14,'FontWeight','bold')
ylabel('Critical gust alleviation factor, K_g','FontSize',14,'FontWeight','bold')
legend('1-DoF aerofoil','FAR25 formula')


