%
% flutterpk.m: Solves dynamic stability of a 2-DoF airfoil using the pk
%              method.
%
% Copyright, Rafael Palacios, July 2022
%            r.palacios@imperial.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all

% Options:
omega_ratio=0.05:0.05:2;    % Ratio of pitch/plunge frequencies.
xalpha=[0 0.05 0.1 0.2];    % Nondimensional position (in c/2) of cg from ea.
Ndof=2;                     % Number of degrees of freedom.
Flag_plot=0;                % =1: Plot root locus for each parameter values.
Tol=1e-3;                   % Convergence tolerance.

% Constants:
rho=1; c=1; omega_a=1;      % Non-dim results are independent of this choice.
mu=5;
r_a =0.25*c;
x_ac=0.25*c;
x_ea=0.35*c;

% Define Theodorsen function.
theod=@(xx) besselh(1,2,xx)./(besselh(1,2,xx)+i*besselh(0,2,xx));


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

    % Obtain GAFs.
    A_0 = [[0 CLa0qs];            [0 0]];
    A_1 = [[CLa0qs CLa1qs+CLa1nc];[0 CMa1qs+CMa1nc]];
    A_2 = [[CLa1nc CLa2nc];       [CMa1nc CMa2nc]];
    A_4 = [[CLa0qs CLa1qs];       [0 0]];
   
    kappa_g=2*[[-1 0];[(x_ea-x_ac)/(c/2) 2]];

    jGAF=0; GAFk=[0:0.1:5]+.0001;  % GAFs at k=0 give problems.
    for kr=GAFk
        jGAF=jGAF+1;
        GAF(jGAF,:,:)=kappa_g*(theod(kr) * (A_0+i*kr*A_4)+ ...
                               i*kr*(A_1-A_4) - kr^2*A_2);
    end

    % mass matrix, and stiffness matrix:
    M=[[m m*(x_cg-x_ea)/(c/2)];[m*(x_cg-x_ea)/(c/2) I_a/(c^2/4)]];
    K=[[k_h 0];[0 k_a/(c^2/4)]];
    

    %% Flutter solver
    Vinf=0.1;  Vmax=1.77; Vdelta=0.01;  % Range of velocities to explore.
    kappa=sqrt(eig(K,M))*c/(2*Vinf);    % Initial guess of the (nondim) eigs.
    
    % Loop through the airspeeds up the diverspeed speed.
    kv=0; Vflutter(kxalpha,komega)=1.78;
    while (Vinf<=Vmax) && (Vflutter(kxalpha,komega)==1.78)
       kv=kv+1;  

       % Loop through the dofs.
       for ke=1:Ndof
           NotConverged=1;
           while NotConverged

               % Interpolate GAF at the current 
               for i1=1:Ndof
                 for i2=1:Ndof
                   GAFint(i1,i2)=interp1(GAFk,squeeze(GAF(:,i1,i2)), ...
                                         kappa(ke),'makima','extrap');
                 end
               end

               % Compute the sqrt of the eigenvalues and enforce real imag
               % part and sort them by their imaginary part.
               p=sqrt(eig(GAFint*(rho*c^2)/8 - K*(c/(2*Vinf))^2,M));
               p=p.*sign(imag(p));
               [psort,isort]=sort(imag(p));
               p=p(isort);

               % Fix-point on the current eigenvalue.
               if (abs(imag(p(ke))-kappa(ke)) < Tol)
                   NotConverged=0;
               end
               kappa(ke)=imag(p(ke));
           end

           % After convergence, compute the rate of decay and check for
           % stability.
           gamma(ke,1)= imag(p(ke))/real(p(ke));
           if gamma(ke)>0
               Vflutter(kxalpha,komega)=Vinf;
           end
       end

       % Compute dimensional eigenvalues and plot them.
       lambda=((gamma+i).*kappa)*(2*Vinf/c);
        
       if Flag_plot
         figure(komega)
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
       end
       data.V(kv)       =Vinf;
       data.lambda(kv,:)=lambda;   
       Vinf=Vinf+Vdelta;
    end

    % Plot eigenvalues in a V-g plot.
    if Flag_plot
        figure(length(omega_ratio)+komega)
        subplot(2,1,1)
            ylabel('real(2\lambdaV\infty/c)')
            grid on, hold on
        subplot(2,1,2)
            xlabel('V_\infty/(c\omega_\alpha)')
            ylabel('imag(2\lambdaV\infty/c)')
            grid on, hold on
        for k=1:2
            subplot(2,1,1), plot(data.V/(c*omega_a),real(data.lambda(:,k)),'b.','MarkerSize',4)
            subplot(2,1,2), plot(data.V/(c*omega_a),imag(data.lambda(:,k)),'b.','MarkerSize',4)
        end
    end
end
end


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
% eof