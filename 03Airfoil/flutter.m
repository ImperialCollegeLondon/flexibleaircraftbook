%
% flutter.m: Solves dynamic stability of a 2-DoF airfoil in state-space
%            description.
% Copyright, Rafael Palacios, June 2018
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Options:
omega_ratio=0.005:0.005:2;  % Ratio of pitch/plunge frequencies.
xalpha=[0 0.05 0.1 0.2];    % Nondimensional position of cg.

Na=4;                       % Order of the RFA (Na<6)
Flag_plot=0;                % Plot root locus for each parameter values.

% Constants:
rho=1; c=1; omega_a=1;   % Non-dim results are independent of this choice.
mu=5;
r_a =0.25*c;
x_ac=0.25*c;
x_ea=0.35*c;


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
    A_3 = [[CLa0qs CLa1qs];[0 0]];
    
    % State-space matrices for unsteady aero.
    for j=1:Na
      A_a((j-1)*2+1:j*2,(j-1)*2+1:j*2)=-eye(2)*b(j);
      B_a((j-1)*2+1:j*2,1:2)=eye(2);
      C_a(1:2,(j-1)*2+1:j*2)= a(j)*(b(j)*A_3-A_0);
    end
      
   
    % Matrix of geometric coefficients for the generalized forces, 
    % mass matrix, and stiffness matrix:
    kappa_g=2*[[-1 0];[(x_ea-x_ac)/(c/2) 2]];
    M=[[m m*(x_cg-x_ea)/(c/2)];[m*(x_cg-x_ea)/(c/2) I_a/(c^2/4)]];
    K=[[k_h 0];[0 k_a/(c^2/4)]];
    
    
    % Loop through the velocities.
    k=0; Vflutter(kxalpha,komega)=0;
    for Vinf=0.01:0.005:2   
       k=k+1;  
       qinf=1/2*rho*Vinf^2;

       M_ae=(4*Vinf^2/c^2)*M - qinf *kappa_g*A_2;
       C_ae=                 - qinf *kappa_g*(A_1-(1/2)*A_3);
       K_ae=               K - qinf *kappa_g*A_0;

       A_ae=[[zeros(2)        eye(2)     zeros(2,2*Na)];
             [ -M_ae\K_ae    -M_ae\C_ae  M_ae\(qinf*kappa_g*C_a)];
             [zeros(2*Na,2)    B_a          A_a]];
       
       % Compute dimensional eigenvalues:
       lambda=eig(A_ae)*(2*Vinf/c);
              
       if Flag_plot
         figure(1)
         if k==1
           plot(real(lambda),imag(lambda),'rs','MarkerSize',8), hold on
           grid on
           xlabel('real(2\lambdaV\infty/c)','FontSize',14,'FontWeight','bold')
           ylabel('imag(2\lambdaV\infty/c)','FontSize',14,'FontWeight','bold')
         elseif mod(k,10)==0
           plot(real(lambda),imag(lambda),'o')
         else            
           plot(real(lambda),imag(lambda),'.')
         end
         data.V(k)       =Vinf;
         data.lambda(k,:)=lambda;   
       end
       
       % Check if there is a positive root.
       if Vflutter(kxalpha,komega)==0 && max(real(lambda))>0
           Vflutter(kxalpha,komega)=Vinf;
       end
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