%
% turb2dof.m: Continuous turbulence response of 2-DoF airfoil using 
%             Jones's approximation to Sears' function.
%
% Dependencies:
%    theodorsen_ajbj.m: aj & bj coefficients of RFA to Theodorsen.
%
% Copyright, Rafael Palacios, April 2023
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Problem parameters:
rho=1; c=1; omega_a=1;   % Non-dim results are independent of this choice.
mu=5;
r_a =0.25*c;
x_ac=0.25*c;
x_ea=0.35*c;
omega_ratio=0.5;      % Ratio of pitch/plunge frequencies.
xalpha=0.2;           % Nondimensional position of cg (in semichords from e.a.)

lturb=1*c;             % Turbulence length scale.

Na=4;                 % Order of the RFA
Ng=5;

MaxFreq=5;            % Max value of omega/omega_\alpha
DeltaFreq=0.002;      % Delta of omega/omega_\alpha


% Obtain a rational-function approximation to Theodorsen:
[a,b]=theodorsen_ajbj(Na);

% Use Jones' approximation to Sears's function.
[a_g, b_g]=sears_ajbj(Ng);
%a_g=[0.5 0.5]'; b_g=[0.13 1.0]'; Ng=2 % Uncomment for Jones coefficients.

% Von Karman turbulence PSD.
PSDKarman=@(xx) ((2*lturb)/pi)*(1+(8/3)*(2*1.339*lturb*xx).^2) ./ ...
                ((1+ (2*1.339*lturb*xx).^2).^(11/6));


%% Construct system matrices.

% Derived parameters:
x_cg = x_ea+xalpha*c/2;
m    = mu*pi*rho*c^2;
I_a  = m*r_a^2;
k_a  = I_a*omega_a^2;
omega_h=omega_ratio*omega_a;
k_h  = m * omega_h^2;
nu_ea=(x_ea-c/2)/(c/2);

% Aerodynamics influence coefficient matrices:
CLa0qs=2*pi;
CLa1qs=2*pi*(1/2-nu_ea);
CMa1qs=-pi/4;
CLa1nc=pi;
CLa2nc=-pi*nu_ea;
CMa1nc=-pi/4;
CMa2nc=-pi/4*(1/4-nu_ea);

A_0 = [[0 CLa0qs];[0 0]];
A_1 = [[CLa0qs CLa1qs+CLa1nc];[0 CMa1qs+CMa1nc]];
A_2 = [[CLa1nc CLa2nc];[CMa1nc CMa2nc]];
A_3 = [[CLa0qs CLa1qs];[0 0]];

% State-space matrices for unsteady aero.
for j=1:Na
  A_a((j-1)*2+1:j*2,(j-1)*2+1:j*2)=-eye(2)*b(j);
  B_a((j-1)*2+1:j*2,1:2)=eye(2);
  C_a(1:2,(j-1)*2+1:j*2)= a(j)*(b(j)*A_3-A_0);
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

% Loop through the velocities.
k=0;
linestyles={'b:','b--','b-.','b-'};
freqs=[0:DeltaFreq:MaxFreq]; % Physical frequencies for plots.

Vinflist=0.2:0.2:0.8; % Vinflist=0.01:0.01:1.04;
for Vinf=Vinflist
   k=k+1;
   if k>4, linestyles{k}='k-'; end
   
   % System matrices for the current dynamic pressure.
   qinf=1/2*rho*Vinf^2;

   M_ae=(4*Vinf^2/c^2)*M - qinf *kappa_g*A_2;
   C_ae=                 - qinf *kappa_g*(A_1-(1/2)*A_3);
   K_ae=               K - qinf *kappa_g*A_0;

   A_ae=[[zeros(2)        eye(2)     zeros(2,2*Na)              zeros(2,Ng)          ];
         [ -M_ae\K_ae    -M_ae\C_ae M_ae\(qinf*kappa_g*C_a)  M_ae\(qinf*kappa_g*C_g) ];
         [zeros(2*Na,2)    B_a          A_a                       zeros(2*Na,Ng)     ];
         [zeros(Ng,2)      zeros(Ng,2)  zeros(Ng,2*Na)              A_g              ]];

   B_ae=[zeros(2,1); zeros(2,1); zeros(2*Na,1); B_g];
   C_ae=[eye(2) zeros(2,2) zeros(2,2*Na) zeros(2,Ng)];
   
   % Construct the state-space description with no feed-through.
   ss_ae=ss(A_ae,B_ae,C_ae,0);
   
   % Compute bode plot to gust inputs and plot output on dimensional
   % frequencies.
   H=freqresp(ss_ae,freqs*c/(2*Vinf));

   figure(1)
   subplot(2,2,1)
      plot(freqs,squeeze(abs(H(1,1,:))),linestyles{k},'LineWidth',1)
      axis([0 2 0 3]), grid on
      title('w_g/V_\infty to \eta')
      ylabel('magnitude')
      hold on
   subplot(2,2,3)
      plot(freqs,squeeze(angle(H(1,1,:)))*180/pi,linestyles{k},'LineWidth',1)
      axis([0 2 -180 180]), grid on
      ylabel('phase (deg)')
      yticks([-180:90:180])
      xlabel('\omega/\omega_\alpha','FontSize',12,'FontWeight','Bold')
      hold on
   subplot(2,2,2)
      plot(freqs,squeeze(abs(H(2,1,:))),linestyles{k},'LineWidth',1)
      axis([0 2 0 3]), grid on
      title('w_g/V_\infty to \alpha','FontSize',12,'FontWeight','Bold')
      ylabel('magnitude','FontSize',12,'FontWeight','Bold')
      hold on
   subplot(2,2,4)
      plot(freqs,squeeze(angle(H(2,1,:)))*180/pi,linestyles{k},'LineWidth',1)      
      axis([0 2 -180 180]), grid on
      hold on
      ylabel('phase (deg)','FontSize',12,'FontWeight','Bold')
      yticks([-180:90:180])
      xlabel('\omega/\omega_\alpha','FontSize',12,'FontWeight','Bold')
      

   % Compute spectral density in dimensional frequency.
   for j1=1:2
       for j2=1:2
           for j3=1:length(freqs)
             psd_ae(j1,j2,j3)=H(j1,1,j3)*PSDKarman(freqs(j3)*c/(2*Vinf))...
                             *conj(H(j2,1,j3))*c/(2*Vinf);
           end
       end
   end
   
   PSDeta=squeeze(abs(psd_ae(1,1,:)));
   PSDalpha=squeeze(abs(psd_ae(2,2,:)));
   
   % Plot PSD.
   figure(2)
     subplot(3,1,1)
       plot(freqs,PSDKarman(freqs*c/(2*Vinf))*c/(2*Vinf),linestyles{k},'LineWidth',1)
       hold on
       ylabel('\Phi_{w/V}','FontSize',14,'FontWeight','Bold')
       axis([0 2 0 2])
     subplot(3,1,2)
       plot (freqs,PSDeta,linestyles{k},'LineWidth',1)
       hold on
       ylabel('\Phi_\eta','FontSize',14,'FontWeight','Bold')
       axis([0 2 0 3])
     subplot(3,1,3)
       plot (freqs,PSDalpha,linestyles{k},'LineWidth',1)
       hold on       
       xlabel('\omega/\omega_\alpha','FontSize',14,'FontWeight','Bold')
       ylabel('\Phi_\alpha','FontSize',14,'FontWeight','Bold')
       axis([0 2 0 3])
       
    % Show that there is perfect coherence.
    figure(3)
     subplot(1,2,1)
       plot (freqs,squeeze((abs(psd_ae(1,2,:)))).^2./(PSDeta.*PSDalpha),linestyles{k},'LineWidth',1)
       hold on
     subplot(1,2,2)
       plot (freqs,squeeze((abs(psd_ae(2,1,:)))).^2./(PSDeta.*PSDalpha),linestyles{k},'LineWidth',1)
       hold on


    % Compute rms and correlation coefficient.
    sigmaw(k)=sqrt( trapz(freqs,PSDKarman(freqs*c/(2*Vinf))*c/(2*Vinf)) );
    sigma1(k)=sqrt(trapz(freqs,PSDeta));
    sigma2(k)=sqrt(trapz(freqs,PSDalpha));
    sigma12(k)=trapz(freqs,0.5*squeeze(psd_ae(1,2,:))+0.5*squeeze(psd_ae(2,1,:)));
end

% Add legends to the figures.
for k=1:3
    figure(k)
    legend('V_\infty=0.2c\omega_\alpha',...
           'V_\infty=0.4c\omega_\alpha',...
           'V_\infty=0.6c\omega_\alpha',...
           'V_\infty=0.8c\omega_\alpha')
end

% Plot PSDs
% (for a high density plot, use Vinflist=0.01:0.01:1.04 above)
figure(4)
 subplot(2,1,1)
 plot(Vinflist,sigma1,'b-','LineWidth',2)
 hold on
 plot(Vinflist,sigma2,'b--','LineWidth',2)
 ylabel('Output RMS')
 axis([0 1.1 0 4])
 plot([1.04 1.04],[0 100],'k:','LineWidth',1)
 legend('\sigma_\eta','\sigma_\alpha')
 subplot(2,1,2)
 plot(Vinflist,real(sigma12./(sigma1.*sigma2)),'b-','LineWidth',2)
 hold on
 plot([1.04 1.04],[-1 1],'k:','LineWidth',1)
 xlabel('V_\infty/(c\omega_\alpha)')
 ylabel('\rho')
 axis([0 1.1 -1 1])

 
 % Plot equal probability ellipses for the very last case
 s1=sigma1(end);
 s2=sigma2(end);
 r =sigma12(end)/(s1*s2);
 

 % eof