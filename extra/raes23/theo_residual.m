% theo_residual.m
%
%  Compute Loewner approximation with polynomial preconditioning of the
%  Theodorsen problem with 3 degrees of freedom.
%
%  It uses the model in Chapter 3 of Palacios & Cesnik (CUP, 2023)
%   https://doi.org/10.1017/9781108354868
%
%
%  This results appeared in a RAeS'23 conference paper
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: Sept 2023. Needs Matlab2023b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
addpath('../../03Airfoil')

% Compute GAFs for 2-DoF airfoil with a flap.
c=1;                  % Not needed, but to make it clearer!
x_ea=0.35*c;          % Elastic axis, from the LE
x_fh=0.85*c;          % Flap hinge location, from the LE
nu_ea=(x_ea-c/2)/(c/2);
nu_fh=(x_fh-c/2)/(c/2);
theta_fh=acos(-nu_fh);

% Aerodynamics influence coefficient matrices 
% (Eq. 3.37 & 3.45 from Palacios & Cesnik, 2023)
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
A_3 = [[0 CLa1qs CLd1qs];
       [0      0      0]];
A_4 = [[CLa0qs CLa1qs CLd1qs];
       [0      0      0]];


%% Convergence of polynomial preconditioning matrices.
if 0
    % Plot k_max dependency.
    j1=0;
    krel=1e-3;
    for klarge=1:0.1:100
        ksmall=klarge*krel;
        clear A H
        j1=j1+1;
    
        k=[klarge-ksmall klarge klarge+ksmall];
        for j=1:length(k)
            ik=1i*k(j);
            H(j,:,:)=(theodorsen(k(j))-1)*A_3 ...
                    +(theodorsen(k(j))-1/2)*ik*A_4;
            A(j,:,:)=A_0+ik*A_1-1/2*ik*A_4+ik^2*A_2 ...
                    +(theodorsen(k(j))-1)*A_3 ...
                    +(theodorsen(k(j))-1/2)*ik*A_4;
        end
        A_0i=reshape(A(1,:,:),2,3);
        
        A_n2=reshape(A(end-2,:,:),2,3); k_n2=klarge-ksmall;
        A_n1=reshape(A(end-1,:,:),2,3); k_n1=klarge;
        A_n0=reshape(A(end,:,:),2,3);   k_n0=klarge+ksmall;
        A_2i=real(-(1/2)*(A_n0-2*A_n1+A_n2)/(ksmall^2));
        A_2ik=-(1/2)*(A_n0-2*A_n1+A_n2)/(ksmall^2);

        A_1bari= real((A_n0-A_n2)/(2*1i*ksmall)) ...
               - 0.5*real((k_n0^2-k_n2^2)*A_2ik/(2*1i*ksmall));
        A_1bar=A_1-(1/2)*A_4;
    
        kplot(j1)=klarge;
        A1plot(j1)=max(max(abs(A_1bari-A_1bar)))/max(max(abs(A_1bar)));
        A2plot(j1)=max(max(abs(A_2i-A_2)))/max(max(abs(A_2)));
    
    end
    figure
    loglog(kplot,A1plot,'k-'), hold on
    loglog(kplot,A2plot,'k--')
    legend('error A_1','error A_2')
    xlabel('k_{max}')
    ylabel('relative error')

end



%% Build a sample of matrices
ksmall=1e-3;
klarge=100;    % Very large to reduce asymptotic errors, but this is 
               % because we are comparing with the analytical solution.
kres=5;
k_n2=klarge+ksmall;
k_n1=klarge+2*ksmall;
k_n0=klarge+3*ksmall;

k=[ksmall 0.05:0.05:kres k_n2 k_n1 k_n0];
for j=1:length(k)
    ik=1i*k(j);
    A(j,:,:)=A_0+ik*A_1-1/2*ik*A_4+ik^2*A_2 ...
            +(theodorsen(k(j))-1)*A_3 ...
            +(theodorsen(k(j))-1/2)*ik*A_4;
end

% Compute polynomial preconditioning matrices.
A_0i=reshape(A(1,:,:),2,3);

A_n2=reshape(A(end-2,:,:),2,3); 
A_n1=reshape(A(end-1,:,:),2,3); 
A_n0=reshape(A(end,:,:),2,3);   
A_2i=real(-(1/2)*(A_n0-2*A_n1+A_n2)/(ksmall^2));

A_2ik=-(1/2)*(A_n0-2*A_n1+A_n2)/(ksmall^2);
A_1bari= real((A_n0-A_n2)/(2*1i*ksmall)) ...
       - 0.5*real((k_n0^2-k_n2^2)*A_2ik/(2*1i*ksmall));

% Obtain samples of approximate residual TF.
numk=length(k)-2;
if mod(numk,2) == 1
  error('need even number of modes')
  sdf
end
for j=1:numk
    ik=1i*k(j);
    H(j,:,:)=(theodorsen(k(j))-1)*A_3 ...
            +(theodorsen(k(j))-1/2)*ik*A_4;
    Hi_temp=reshape(A(j,:,:),2,3)-A_0i-ik*A_1bari-(ik^2)*A_2i;
    Hi(j,1:2,1:3)=Hi_temp(1:2,1:3);
    Hi(j,3,1:3)=Hi_temp(2,1:3); % Added to have a rectangular matrix.
end

% The very last matrix is the feedthrough gain. Remove
Hi_inf(1:2,1:3)=Hi_temp(1:2,1:3);
Hi_inf(3,1:3)=Hi_temp(2,1:3); % Added to have a rectangular matrix.
%Hi_inf=(1/8)*A_4-(1/2)*A_3;
for j=1:numk
    Hi(j,:,:)=reshape(Hi(j,:,:),3,3)-Hi_inf;
end


% Plot the approximated matrix.
if 0
    maxkplot=20;
    figure
    subplot(2,2,1)
      plot(k(1:numk),squeeze(abs(Hi(:,1,1))),'k','LineWidth',2)
      title('H_{11}')
      xlim([0 maxkplot])
      ylabel('Magnitude (abs)')
    subplot(2,2,3)
       plot(k(1:numk),squeeze(angle(Hi(:,1,1)))*180/pi,'k','LineWidth',2)
       xlim([0 maxkplot])
       ylim([0 180]), yticks([0 90 180]) 
       ylabel('Phase (deg)'), xlabel('k')
    subplot(2,2,2)
      plot(k(1:numk),squeeze(abs(Hi(:,1,2))),'k','LineWidth',2)
      xlim([0 maxkplot])
      title('H_{12}')
    subplot(2,2,4)
      plot(k(1:numk),squeeze(angle(Hi(:,1,2)))*180/pi,'k','LineWidth',2)
      ylim([-180 0]), yticks([-180 -90 0])
      xlim([0 maxkplot]), xlabel('k')
end 

% Loewner matrix.
%Split between left and right data.
nrows=3;
ncols=nrows; % Need a bit of work to take separate numbers. For this example
             % we just repeat one row.
Vl=[]; Wr=[];
for ii=1:numk
    if mod(ii,2)==1 % odd numbers
      Ml(ii)  = 1i*k(ii);
      Ml(ii+1)=-1i*k(ii);
      Vl(end+1:end+nrows,:)= (reshape(Hi(ii,:,:),nrows,ncols));
      Vl(end+1:end+nrows,:)= (reshape(conj(Hi(ii,:,:)),nrows,ncols));
    else
      Mr(ii-1)= 1i*k(ii);
      Mr(ii)  =-1i*k(ii);
      Wr(:,end+1:end+ncols)= (reshape(Hi(ii,:,:),nrows,ncols));
      Wr(:,end+1:end+ncols)= (reshape(conj(Hi(ii,:,:)),nrows,ncols));
    end
end

% Construct Loewner matrix
 for ii=1:numk
     i1=nrows*(ii-1);
     for jj=1:numk
         j1=ncols*(jj-1);
         Loew(i1+1:i1+nrows,j1+1:j1+ncols)= (1/(Ml(ii)-Mr(jj))) * ...
             ((Vl(i1+1:i1+nrows,:)) - Wr(:,j1+1:j1+ncols));
         Loes(i1+1:i1+nrows,j1+1:j1+ncols)= (1/(Ml(ii)-Mr(jj))) * ...
             (Ml(ii)*(Vl(i1+1:i1+nrows,:)) ...
             -Mr(jj)*Wr(:,j1+1:j1+ncols));
     end
 end

% Tranform into a real system.
s2=1/sqrt(2);
for ii=1:numk/2
    i1=(2*nrows)*(ii-1);
    Jl(i1+1:i1+nrows,i1+1:i1+nrows)=s2*eye(nrows);
    Jl(i1+1:i1+nrows,i1+nrows+1:i1+2*nrows)=-1i*s2*eye(nrows);
    Jl(i1+nrows+1:i1+2*nrows,i1+1:i1+nrows)= s2*eye(nrows);
    Jl(i1+nrows+1:i1+2*nrows,i1+nrows+1:i1+2*nrows)= 1i*s2*eye(nrows);
end

Loew_re=real(Jl'*Loew*Jl);
Loes_re=real(Jl'*Loes*Jl);
Vl_re=real(Jl'*Vl);
Wr_re=real(Wr*Jl);

% SVD projection
[U1,S1,V1]=svd([Loew_re Loes_re]);
[U2,S2,V2]=svd([Loew_re;Loes_re]);
Nmax=5;  % Defines the size of the ROM.
E_loew=-transpose(U1(:,1:Nmax))*Loew_re*V2(:,1:Nmax);
A_loew=-transpose(U1(:,1:Nmax))*Loes_re*V2(:,1:Nmax);
B_loew=transpose(U1(:,1:Nmax))*Vl_re;
C_loew=Wr_re*V2(:,1:Nmax);
D_loew=real(Hi_inf);
sys_loew=dss(A_loew,B_loew,C_loew,D_loew,E_loew);

% Plot results.
figure
kappai=0.01:0.01:kres/2;
ii=0;
for i1=1:2
  for i2=1:3
    ii=ii+1;
    GAF_loew=squeeze(freqresp(sys_loew(i1,i2),kappai));
    subplot(2,3,ii)
        plot(k(1:5:numk),real(squeeze(H(1:5:end,i1,i2))),'sb'), hold on
        plot(k(1:5:numk),imag(squeeze(H(1:5:end,i1,i2))),'ob'), hold on
        plot(kappai,real(GAF_loew),'k-')
        plot(kappai,imag(GAF_loew),'k--')
        xlim([0 2])
  end
end
subplot(2,3,1)
  ylabel('H_{L\eta}','FontSize',12,'FontWeight','bold')
subplot(2,3,2)
  ylabel('H_{L\alpha}','FontSize',12,'FontWeight','bold')
subplot(2,3,3)
  ylabel('H_{L\delta}','FontSize',12,'FontWeight','bold')
  legend('Exact, real', 'Exact, imag','Approx, real','Approx, imag')
subplot(2,3,4)
  ylabel('H_{M\eta}','FontSize',12,'FontWeight','bold')
  xlabel('k','FontSize',12,'FontWeight','bold')
subplot(2,3,5)
  ylabel('H_{M\alpha}','FontSize',12,'FontWeight','bold')
  xlabel('k','FontSize',12,'FontWeight','bold')
subplot(2,3,6)
  ylabel('H_{M\delta}','FontSize',12,'FontWeight','bold')
  xlabel('k','FontSize',12,'FontWeight','bold')

% Check eigs.
if (max(real(eig(sys_loew))) > 0)
    error('Unstable system')
end
% eof