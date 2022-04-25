%
% turb2dof_lqq.m: Continuous turbulence response of 2-DoF airfoil in 
%                 closed-loop with LQG control.
%
% Copyright, Rafael Palacios, October 2020
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Problem parameters:
rho=1; c=1; omega_a=1;    % Results are independent of this choice. This
                          % is kept to better understand the problem.
mu=5;
r_a =0.25*c;
x_ac=0.25*c;
x_ea=0.35*c;
x_fh=0.85*c;          % Flap hinge location.
omega_ratio=0.5;      % Ratio of plunge/pitch frequencies.
xalpha=0.2;           % Nondimensional position of cg from e.a.
lturb=1*c;            % Turbulence length scale.
aturb=1.339;          % Von Karman constant.
Na=4;                 % Order of the RFA (Na<6)

Vinf=0.5*c*omega_a;
qinf=1/2*rho*Vinf^2;
sigma_w=0.1*Vinf;     % gust intensity.

Qmax =100;             % Control weights.

freqs=[0.01:.002:2]*omega_a; % Physical frequencies for plots.


%% Obtain a rational-function approximation to Theodorsen:
theod=@(xx) besselh(1,2,xx)./(besselh(1,2,xx)+i*besselh(0,2,xx));
deltak=0.01;
k=0:deltak:5;
Wk=ones(size(k));
Wk(1:1/deltak)=100;
systh=frd(theod(k)-0.5,k);
systh4=fitmagfrd(systh,Na,1,Wk)+0.5;

% Identify Jones-type approximation.
b=-pole(systh4);
pol=poly(pole(systh4))-poly(zero(systh4))/2;
 
 A=zeros(Na);
 for i1=1:Na
     A(1,i1)=A(1,i1)+1;
     for i2=1:Na, if i2 ~= i1
      A(2,i1)=A(2,i1)+b(i2);
     for i3=1:Na, if i3<i2 && i3~=i1
      A(3,i1)=A(3,i1)+b(i2)*b(i3);
     for i4=1:Na, if i4<i3 && i4<i2 && i4~=i1
      A(4,i1)=A(4,i1)+b(i2)*b(i3)*b(i4);
     for i5=1:Na, if i5<i4 && i5<i3 && i5<i2 && i5~=i1
      A(5,i1)=A(5,i1)+b(i2)*b(i3)*b(i4)*b(i5);
     end, end, end, end, end, end, end, end
 end
a=A\pol(1:Na)';
% Check that the sum of a is 0.5, which is not enforced.
if abs(sum(a)-0.5) > 0.01
    stop 'Needs a better RFA'
end
clear A pol

%% Use Jones' approximation to Sears.
a_g=[0.5 0.5]';
b_g=[0.13 1.0]';
Ng=2;

%% Generate  Von Karman vertical turbulence PSD (unit T.I.) - Eq. (2.44)
% Write in terms of reduced frequency.
PSDKarman=@(xx) (lturb/pi)*(1+(8/3)*((2*aturb*lturb/c)*xx).^2) ...
                       ./ ((1+      ((2*aturb*lturb/c)*xx).^2).^(11/6));

                   
%% Construct system matrices.
% Derived parameters:
x_cg = x_ea+xalpha*c/2;   % CG coordinate
m    = mu*pi*rho*c^2;     % Physical mass
I_a  = m*r_a^2;           % Moment of inertia about e.a.
k_a  = I_a*omega_a^2;
omega_h=omega_ratio*omega_a;
k_h  = m * omega_h^2;
nu_ea=(x_ea-c/2)/(c/2);
nu_fh=(x_fh-c/2)/(c/2);
theta_fh=acos(-nu_fh);

% Aerodynamics influence coefficient matrices 
% (Eq. 3.37 & 3.45)
CLa0qs=2*pi;
CLa1qs=2*pi*(1/2-nu_ea);
CMa1qs=-pi/4;
CLa1nc=pi;
CLa2nc=-pi*nu_ea;
CMa1nc=-pi/4;
CMa2nc=-pi/4*(1/4-nu_ea);

% (Eq. 3.38 & 3.46)
CLd0qs= 2*pi - 2*theta_fh + 2*sin(theta_fh);
CLd1qs= (1/2-nu_fh)*(2*pi-2*theta_fh) ...
      + (2-nu_fh)*sin(theta_fh);
CMd0qs=-(1/2)*(1+nu_fh)*sin(theta_fh);
CMd1qs=-(1/4)*(pi-nu_fh + (2/3)*(1/2-nu_fh)*(2+nu_fh)*sin(theta_fh));
CLd1nc=pi-theta_fh-nu_fh*sin(theta_fh);
CLd2nc=-nu_fh*(pi-theta_fh) + (1/3)*(2+nu_fh^2)*sin(theta_fh);
CMd1nc=-(1/4)*(pi-theta_fh+(2/3-nu_fh-(2/3)*nu_fh^2)*sin(theta_fh));
CMd2nc=-(1/4)*((1/4-nu_fh)*(pi-theta_fh) ...
              +(2/3-(5/12)*nu_fh+(1/3)*nu_fh^2+(1/6)*nu_fh^3) ...
               *sin(theta_fh));

% Eq. (3.60)          
A_0 = [[0 CLa0qs CLd0qs];
       [0 0      CMd0qs]];
A_1 = [[CLa0qs CLa1qs+CLa1nc CLd1qs+CLd1nc];
       [0      CMa1qs+CMa1nc CMd1qs+CMd1nc]];
A_2 = [[CLa1nc CLa2nc CLd2nc];
       [CMa1nc CMa2nc CMd2nc]];
A_3 = [[CLa0qs CLa1qs CLd1qs];
       [0      0      0]];

% State-space matrices for unsteady aero - Eq. (3.69)
for j=1:Na
  A_a((j-1)*3+1:j*3,(j-1)*3+1:j*3)=-eye(3)*b(j);
  B_a((j-1)*3+1:j*3,1:3)=eye(3);
  C_a(1:2,(j-1)*3+1:j*3)= a(j)*(b(j)*A_3-A_0);
end

% State-space matrices for gust aero.
for j=1:Ng
   A_g(j,j)=-b_g(j);
   B_g(j,1)=1;
   C_g(1:2,j)=a_g(j)*b_g(j)*A_0(1:2,2);
end

% Matrix of geometric coefficients for the generalized forces, 
% mass matrix, and stiffness matrix:
kappa_g=2*[[-1 0];[(x_ea-x_ac)/(c/2) 2]];
M=[[m m*(x_cg-x_ea)/(c/2)];[m*(x_cg-x_ea)/(c/2) I_a/(c^2/4)]];
K=[[k_h 0];[0 k_a/(c^2/4)]];


% Obtain time history of Gaussian noise. Use same one for all controllers.
Ntimes=4000;               % Number of time steps.
deltaw=randn(Ntimes+1,1);  % White noise.
dt=0.05; t=0:dt:Ntimes*dt; % Non-dimensional times for simulation.


%% Loop through the Q matrices
k=0;
linestyles={'b--','b:','b-'};
markstyles={'bo','bx','bs'};

for k=1:3
   % Define LQR cost.
   Nx=3*(Na+2)+Ng+4;
   if k>1
       Q=diag(ones(Nx,1));
       switch k
           case 2, Q(1,1)=Qmax;
           case 3, Q(2,2)=Qmax;
       end
   end
  
   %% Aeroelastic system
   % System matrices for the current dynamic pressure.
   % Equations are written in reduced frequencies.

   M_ae=(4*Vinf^2/c^2)*M - qinf *kappa_g*A_2(:,1:2);
   invM=inv(M_ae);
   C_ae= - qinf *kappa_g*(A_1-(1/2)*A_3);
   K_ae= - qinf *kappa_g*A_0;
   K_ae(1:2,1:2)=K_ae(1:2,1:2)+ K;

   % DoFs in A_ae: eta, alpha, delta, eta', alpha', delta', x_a, x_g
   A_ae=[[zeros(3)        eye(3)      zeros(3,3*Na)          zeros(3,Ng)          ];
         [ -invM*K_ae    -invM*C_ae   invM*qinf*kappa_g*C_a  invM*qinf*kappa_g*C_g];
         [zeros(1,3)      zeros(1,3)  zeros(1,3*Na)          zeros(1,Ng)          ];
         [zeros(3*Na,3)   B_a         A_a                    zeros(3*Na,Ng)       ];
         [zeros(Ng,3)     zeros(Ng,3) zeros(Ng,3*Na)         A_g                  ]];

   % inputs: \ddot{delta}, w_g
   B_ae=[[zeros(3,1)                    zeros(3,1)];
         [invM*qinf*kappa_g*A_2(:,3)    zeros(2,1)];
         [1                             0         ];
         [zeros(3*Na,1)                 zeros(3*Na,1)]; 
         [zeros(Ng,1)                   B_g]];
   C_ae=[[eye(3) zeros(3,3*(Na+2)+Ng-3)];
        [ zeros(1,4) 1 zeros(1,3*(Na+2)+Ng-5)]];
   
   % System with two inputs and four outputs.
   ss_ae=ss(A_ae,B_ae,C_ae,0);
   
   %% Turbulence filter
   % Campbell Rational function approximation, in terms of omega*T.
   T=2*aturb*lturb/c;
   GCampbell=sqrt(T/(pi*aturb)) ...
            *tf([sqrt(8/3)*T 1],[T 1]) ...
            *tf([91/12*(T^2) 52*T 60],[935/216*(T^3) 561/12*(T^2) 102*T 60]);
   ssCampbell=append(1,ss(GCampbell,'minimal'));
   
   % connect turbulent filter to aeroelastic system.
   ss_t=series(ssCampbell,ss_ae);
   
   if k>1
       % Create LQG. Measurement is pitch rate. Output is \ddot{flap}.
       [k1,L1,p1]=kalman (ss_t(4,:),diag([1e-8 1e-2]),1e-6);
       [Klqr,Slqr,Elqr] = lqr(ss_t(:,1),Q,1);
       A_lqg=[ss_t.A-ss_t.B(:,1)*Klqr ss_t.B(:,1)*Klqr;
              zeros(size(ss_t.A)) ss_t.A-L1*ss_t.C(4,:)];
       B_lqg=[ss_t.B(:,2); ss_t.B(:,2)];
       C_lqg=[ss_t.C zeros(size(ss_t.C))];
       
       ss_ae=ss(A_lqg, B_lqg, C_lqg, 0);
   else
       ss_ae=ss_t(:,2);  %Ignore control.
   end
   
   % Plot eigenvalues (given in non-dimensional time -reduced frequency).
   V_ae=eig(ss_ae);
   figure(1)
   plot(real(V_ae),imag(V_ae),markstyles{k},'MarkerSize',10), hold on
   xlabel('real(\lambda)'), grid on
   ylabel('imag(\lambda)')
   axis([-1.2 0 -1.5 1.5])
   
   %% Compute bode plot to gust inputs on in reduced frequencies.
   H_ae=freqresp(ss_ae,freqs*c/(2*Vinf));

   % Plot admittances in terms of the physical frequencies.
   figure(2)
   subplot(2,3,1)
      plot(freqs/omega_a,squeeze(abs(H_ae(1,1,:))),linestyles{k},'LineWidth',1), hold on
      axis([0 1.5 0 1.5]), grid on
      title('w_g/V_\infty to \eta')
      ylabel('magnitude')
   subplot(2,3,4)
      plot(freqs/omega_a,squeeze(angle(H_ae(1,1,:)))*180/pi,linestyles{k},'LineWidth',1), hold on
      axis([0 1.5 -180 180]), grid on
      ylabel('phase (deg)')
      yticks([-180:90:180])
      xlabel('\omega/\omega_\alpha','FontSize',12,'FontWeight','Bold')
   subplot(2,3,2)
      plot(freqs/omega_a,squeeze(abs(H_ae(2,1,:))),linestyles{k},'LineWidth',1), hold on
      axis([0 1.5 0 1.5]), grid on
      title('w_g/V_\infty to \alpha','FontSize',12,'FontWeight','Bold')
      ylabel('magnitude')
   subplot(2,3,5)
      plot(freqs/omega_a,squeeze(angle(H_ae(2,1,:)))*180/pi,linestyles{k},'LineWidth',1), hold on
      axis([0 1.5 -180 180]), grid on     
      ylabel('phase (deg)')
      yticks([-180:90:180])
      xlabel('\omega/\omega_\alpha','FontSize',12,'FontWeight','Bold')
   subplot(2,3,3)
      plot(freqs/omega_a,squeeze(abs(H_ae(3,1,:))),linestyles{k},'LineWidth',1), hold on
      axis([0 1.5 0 1.5]), grid on
      title('w_g/V_\infty to \delta')
      ylabel('magnitude')
   subplot(2,3,6)
      plot(freqs/omega_a,squeeze(angle(H_ae(3,1,:)))*180/pi,linestyles{k},'LineWidth',1), hold on
      axis([0 1.5 -180 180]), grid on
      ylabel('phase (deg)')
      yticks([-180:90:180])
      xlabel('\omega/\omega_\alpha','FontSize',12,'FontWeight','Bold')
      
   %% Compute the spectral density in dimensional angular frequency.
   for j1=1:3
       for j2=1:2
           for j3=1:length(freqs)
             psd(j1,j2,j3)=H_ae(j1,1,j3)*conj(H_ae(j2,1,j3))*c/(2*Vinf);
           end
       end
   end  
   PSDeta=sigma_w*squeeze(abs(psd(1,1,:)))/2;
   PSDalpha=sigma_w*squeeze(abs(psd(2,2,:)))/2;
   
   % Plot PSD.
   figure(3)
     subplot(3,1,1)
       plot(freqs/omega_a,sigma_w*PSDKarman(freqs*c/(2*Vinf))*c/(2*Vinf),...
            linestyles{k},'LineWidth',1)
       hold on
       ylabel('\Phi_{w}','FontSize',14,'FontWeight','Bold')
       axis([0 1.5 0 0.02])
     subplot(3,1,2)
       plot (freqs/omega_a,PSDeta,linestyles{k},'LineWidth',1)
       hold on
       ylabel('\Phi_\eta','FontSize',14,'FontWeight','Bold')
       axis([0 1.5 0 0.03])
     subplot(3,1,3)
       plot (freqs/omega_a,PSDalpha,linestyles{k},'LineWidth',1)
       hold on       
       xlabel('\omega/\omega_\alpha','FontSize',14,'FontWeight','Bold')
       ylabel('\Phi_\alpha','FontSize',14,'FontWeight','Bold')
       axis([0 1.5 0 0.01])
       
%    % Compute rms and correlation coefficient.
%    % It can be checked that trapz{PSDKarmanX+2PSDKarman} gives sigma_w.
%     sigmaw(k)=sqrt( trapz(freqs,PSDKarman(freqs*c/(2*Vinf)))*c/(2*Vinf) );
%     sigma1(k)=sqrt(trapz(freqs,PSDeta));
%     sigma2(k)=sqrt(trapz(freqs,PSDalpha));
%     sigma12(k)=trapz(freqs,0.5*squeeze(psd(1,2,:))+0.5*squeeze(psd(2,1,:)));



    %% Compute time history to gust inputs.
    % Simulate time response of the compound system to the white noise. 
    [q,t] =lsim(ss_ae,deltaw,t,'foh');
    t1=t*omega_a*c/(2*Vinf);
    figure (4)
    subplot(3,1,1)
       plot(t1,q(:,1),linestyles{k},'LineWidth',1)
       hold on
       ylabel('\eta','FontSize',14,'FontWeight','Bold')
    subplot(3,1,2)
       plot(t1,q(:,2)*180/pi,linestyles{k},'LineWidth',1)
       hold on
       ylabel('\alpha','FontSize',14,'FontWeight','Bold')
    subplot(3,1,3)
       plot(t1,q(:,3)*180/pi,linestyles{k},'LineWidth',1)
       hold on
       ylabel('\delta','FontSize',14,'FontWeight','Bold')
       xlabel('\omega_\alphat','FontSize',14,'FontWeight','Bold')
    
end

% Add legends to the figures.
for k=1:4
    figure(k)
    legend('Open loop',...
           ['LQG (Q_{11}=' int2str(Qmax) ')'],...
           ['LQG (Q_{22}=' int2str(Qmax) ')'],'FontSize',12)
end
 
% eof