%
% flutter.m: Solves dynamic stability of a 2-DoF airfoil in state-space
%            description. Example 3.3 of Palacios & Cesnik's book.
%
% Copyright, Rafael Palacios, May 2023
%            r.palacios@imperial.ac.uk
%
% Dependencies:
%    theodorsen_ajbj.m: Coefficients of the RFA of Theodorsen's lift 
%                       deficiency function.
%
% Note: For solution in frequency-domain, see flutterpk.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Parametric variations:
omega_ratio=0.005:0.005:2;  % Ratio of pitch/plunge frequencies.
xalpha=[0 0.05 0.1 0.2];    % Nondimensional position of cg.

% Options.
Na=4;                       % Order of the RFA.
Flag_plot=0;                % Plot root locus for each parameter values.

% Problem constants:
rho=1; c=1; omega_a=1;     % Non-dim results are independent of this choice.
mu=5;
r_a =0.25*c;
x_ac=0.25*c;
x_ea=0.35*c;

% Obtain a rational-function approximation to Theodorsen:
a=zeros(1,Na);
b=zeros(1,Na);
[a,b]=theodorsen_ajbj(Na)


%% Loop through the parameters in the problem.
for komega=1:length(omega_ratio)
for kxalpha=1:length(xalpha)
    omega_ratio(komega)

    % Derived parameters:
    x_cg = x_ea+xalpha(kxalpha)*c/2;
    m    = mu*pi*rho*c^2;
    I_a  = m*r_a^2;
    k_a  = I_a*omega_a^2;
    omega_h=omega_ratio(komega)*omega_a;
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
    A_4 = [[CLa0qs CLa1qs];[0 0]];
    
    % State-space matrices for unsteady aero.
    for j=1:Na
      A_a((j-1)*2+1:j*2,(j-1)*2+1:j*2)=-eye(2)*b(j);
      B_a((j-1)*2+1:j*2,1:2)=eye(2);
      C_a(1:2,(j-1)*2+1:j*2)= a(j)*(b(j)*A_4-A_0);
    end
      
   
    % Matrix of geometric coefficients for the generalized forces, 
    % mass matrix, and stiffness matrix:
    kappa_g=2*[[-1 0];[(x_ea-x_ac)/(c/2) 2]];
    M=[[m m*(x_cg-x_ea)/(c/2)];[m*(x_cg-x_ea)/(c/2) I_a/(c^2/4)]];
    K=[[k_h 0];[0 k_a/(c^2/4)]];
    
    
    % Loop through the vthrough the airspeeds up the diverspeed speed.
    Vinf=0.01;  Vmax=1.77; Vdelta=0.005;  % Range of velocities to explore.
    kv=0; Vflutter(kxalpha,komega)=1.78;
    while (Vinf<=Vmax) && (Vflutter(kxalpha,komega)==1.78)
       kv=kv+1;  
       qinf=1/2*rho*Vinf^2;

       M_ae=(4*Vinf^2/c^2)*M - qinf *kappa_g*A_2;
       C_ae=                 - qinf *kappa_g*(A_1-(1/2)*A_4);
       K_ae=               K - qinf *kappa_g*A_0;

       A_ae=[[zeros(2)        eye(2)     zeros(2,2*Na)];
             [ -M_ae\K_ae    -M_ae\C_ae  M_ae\(qinf*kappa_g*C_a)];
             [zeros(2*Na,2)    B_a          A_a]];
       
       % Compute dimensional eigenvalues:
       lambda=eig(A_ae)*(2*Vinf/c);
              
       if Flag_plot
         figure(1)
         if kv==1
           plot(real(lambda),imag(lambda),'rs','MarkerSize',8), hold on
           grid on
           xlabel('real(2\lambdaV\infty/c)','FontSize',14,'FontWeight','bold')
           ylabel('imag(2\lambdaV\infty/c)','FontSize',14,'FontWeight','bold')
         elseif mod(kv,10)==0
           plot(real(lambda),imag(lambda),'o')
         else            
           plot(real(lambda),imag(lambda),'.')
         end
         data.V(kv)       =Vinf;
         data.lambda(kv,:)=lambda;   
       end
       
       % Check if there is a positive root.
       if max(real(lambda))>0
           Vflutter(kxalpha,komega)=Vinf;
       end
       Vinf=Vinf+Vdelta;
    end

    % V-g plot.
    if Flag_plot
        figure
        subplot(2,1,1)
            ylabel('real(2\lambdaV\infty/c)')
            grid on, hold on
        subplot(2,1,2)
            xlabel('V_\infty/(c\omega_\alpha)')
            ylabel('imag(2\lambdaV\infty/c)')
            grid on, hold on
        for k=1:4+2*Na
            subplot(2,1,1), plot(data.V/(c*omega_a),real(data.lambda(:,k)),'b.','MarkerSize',4)
            subplot(2,1,2), plot(data.V/(c*omega_a),imag(data.lambda(:,k)),'b.','MarkerSize',4)
        end
    end
end
end


if Flag_plot ==0
    % Create flutter envelope plot.
    figure
    plot(omega_ratio,Vflutter(1,:),'b--','LineWidth',3), hold on
    plot(omega_ratio,Vflutter(2,:),'b-','LineWidth',3)
    plot(omega_ratio,Vflutter(3,:),'b:','LineWidth',3)
    plot(omega_ratio,Vflutter(4,:),'b-.','LineWidth',3)
    grid on
    xlabel('\omega_h/\omega_\alpha','FontSize',16,'FontWeight','bold')
    ylabel('V_{flutter}/(c\omega_\alpha)','FontSize',16,'FontWeight','bold')
    axis([0 2 0 2])
    legend('x_{cg}=0.35c','x_{cg}=0.375c','x_{cg}=0.40c','x_{cg}=0.45c')
end
% eof