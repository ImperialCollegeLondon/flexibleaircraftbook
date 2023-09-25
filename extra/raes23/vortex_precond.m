% vortex_precond.m
%
%  Loewner interpolation with frequency-limited preconditioning for the 
%  state-spaceformulation of the 2D airfoil unsteady aerodynamics using the
%  discrete vortex method.
%
%  It has been modified from Example 7.2 in Palacios & Cesnik (CUP, 2023)
%   https://doi.org/10.1017/9781108354868
%
%  A single dependency on each of the 3 degree of freedom is now included, 
%  with derivatives introduce through a second-oder operator.
%
%  This results appeared in a RAeS'23 conference paper
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: Sept 2023. Needs Matlab2023b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all
addpath('../../07UnsteadyAero')

Nb=50;                         % Number of segments on aerofoil
Nw=Nb*30;                      % Number of segments in wake 
dx=1/Nb;                       % Nondimensional panel length
ds=2*dx;                       % Nondimensional time step.
knyq=pi/ds;                    % Nyquist reduced frequency.


% Define non-dimensional coordinates from leading edge.
x=dx/4:dx:1+Nw/Nb-((3*dx)/4);  % Coordinates of vortices (aerofoil/wake).
xi=3*dx/4:dx:1-dx/4;           % Coordinates of the Nb collocation points.
                    
x_ea=0.35;                     % Pitch rotations about 0.35c.
x_fh=0.85;                     % Flap hinge at 0.85c.

linestyles={'--',':','-.'};

%% Aerofoil system equations in discrete time.
[Aa,Ba,Ca,Da]=vortex_getsys(Nb,Nw,x-1/4);
sysa=ss(Aa,Ba,Ca/2/pi,Da/2/pi,ds);

% Contribution of airfoil motions to input matrix.
z=tf('z',ds);
D1=1/(2*ds)*(1/z^2)*(3*z^2-4*z+1);

W(1:Nb,1)= ones(Nb,1)+D1*(2*(xi'-x_ea));
W(1:Nb,2)= D1*ones(Nb,1);
for j=1:Nb
    if xi(j) >x_fh
        W(j,3)=1+D1*(2*(xi(j)-x_fh));
    else
        W(j,3)=0;
    end
end

if 1
    % Build the full system and its rrequency-limited balanced approximation. 
    sys=sysa*W;
    R2 = reducespec(sys,"balanced");
    R2.Options.FreqIntervals = [0,5];
    R2 = process(R2);
    rsys = getrom(R2,Order=8,Method="matchDC");
    save('rsys.mat','rsys')
else
    load rsys
end



%% Compute polynomial preconditioning matrices.
ksmall=1e-3;
klarge=2;
k_n2=klarge-ksmall;
k_n1=klarge;
k_n0=klarge+ksmall;

A_n=freqresp(rsys,[k_n2 k_n1 k_n0]);
A_n2=reshape(A_n(:,:,1),2,3); 
A_n1=reshape(A_n(:,:,2),2,3); 
A_n0=reshape(A_n(:,:,3),2,3);

A_2i=real(-(1/2)*(A_n0-2*A_n1+A_n2)/(ksmall^2));
A_0i=real(reshape(A_n(:,:,2),2,3)+A_2i*klarge^2);
%A_n0=A_n0+k_n0^2*A_2i;
%A_n2=A_n2+k_n2^2*A_2i;
A_1i= real((A_n0-A_n2)/(2*1i*ksmall));


k=[0.001 0.05:0.05:2];
numk=length(k);
A=freqresp(rsys,k);

% Plot polynomial preconditioning vs the actual matrix.
if 1
    figure
    for jj=1:numk
        ik=1i*k(jj);
        Ai(:,:,jj)=A_0i+ik*A_1i+(ik^2)*A_2i;
    end
    figure
    for kk=1:2  % Repeat for each output in the system (lift, moment).
      for jj=1:3 
        subplot(2,3,(kk-1)*3+jj)
           yyaxis left
             plot(k(1:numk),squeeze(real(Ai(kk,jj,1:numk))),'LineWidth',2), hold on
             plot(k(1:numk),squeeze(real(A(kk,jj,1:numk))),':','LineWidth',2), hold on 
             xlabel('k','FontSize',12,'FontWeight','bold')
             ylabel(['A' int2str(kk) int2str(jj)],'FontSize',12,'FontWeight','bold'), hold on
           yyaxis right
             plot(k(1:numk),squeeze(imag(Ai(kk,jj,1:numk))),'LineWidth',2)
             plot(k(1:numk),squeeze(imag(A(kk,jj,1:numk))),':','LineWidth',2), hold on 
      end
    end
end


% Obtain samples of approximate residual TF.
for jj=1:numk
    ik=1i*k(jj);
    Hi_temp=reshape(A(:,:,jj),2,3)-A_0i-ik*A_1i-(ik^2)*A_2i;
    Hi(jj,1:2,1:3)=Hi_temp(1:2,1:3);
    Hi(jj,3,1:3)=Hi_temp(2,1:3); % Added to have a rectangular matrix.
end


% Rescale to one at k=0 and add a zero at a large number
for kk=1:3
   ScalingY(kk,kk)=max(max(abs(Hi(:,kk,:))));
end
for jj=1:3
   ScalingU(jj,jj)=max(max(abs(Hi(:,:,jj))));
end
for jj=1:numk
   Hi(jj,1:3,1:3)=inv(ScalingY)*reshape(Hi(jj,1:3,1:3),3,3)*inv(ScalingU);
end
numk=numk+1;
k=[k 10];
Hi(numk,1:2,1:3)=0;

if 0
    figure
    for kk=1:2  % Repeat for each output in the system (lift, moment).
      for jj=1:3 
        subplot(2,3,(kk-1)*3+jj)
           yyaxis left,  plot(k(1:numk),squeeze(real(Hi(:,kk,jj))),'LineWidth',2)
           xlabel('k','FontSize',12,'FontWeight','bold')
           ylabel(['A' int2str(kk) int2str(jj)],'FontSize',12,'FontWeight','bold')
           yyaxis right, plot(k(1:numk),squeeze(imag(Hi(:,kk,jj))),'LineWidth',2)
      end
    end
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



%% Plot results
figure
% Plot reference.
for kk=1:2  % Repeat for each output in the system (lift, moment).
  for jj=1:3 
    subplot(2,3,(kk-1)*3+jj)
       yyaxis left
         plot(k(1:numk-1),squeeze(real(A(kk,jj,1:numk-1))),'-','LineWidth',2), hold on 
         xlabel('k','FontSize',12,'FontWeight','bold')
         ylabel(['A' int2str(kk) int2str(jj)],'FontSize',12,'FontWeight','bold'), hold on
       yyaxis right
         plot(k(1:numk-1),squeeze(imag(A(kk,jj,1:numk-1))),'-','LineWidth',2), hold on 
  end
end

Nmax=[4 6 8];
kappai=0.01:0.01:2;
iksys=tf([1 0],1);        % Differential operator.
for ii=1:length(Nmax)
    E_loew=-transpose(U1(:,1:Nmax(ii)))*Loew_re*V2(:,1:Nmax(ii));
    A_loew=-transpose(U1(:,1:Nmax(ii)))*Loes_re*V2(:,1:Nmax(ii));
    B_loew=transpose(U1(:,1:Nmax(ii)))*Vl_re;
    C_loew=Wr_re*V2(:,1:Nmax(ii));
    D_loew=[];
    sys_loew=ScalingY*dss(A_loew,B_loew,C_loew,D_loew,E_loew)*ScalingU;
    GAF_loew=squeeze(freqresp(sys_loew(1:2,1:3) ...
                    +A_0i+iksys*A_1i+iksys*iksys*A_2i,kappai));

    for kk=1:2
        for jj=1:3
            subplot(2,3,(kk-1)*3+jj)
            yyaxis left
               plot(kappai,squeeze(real(GAF_loew(kk,jj,:))),linestyles{ii},'LineWidth',2)
            yyaxis right
               plot(kappai,squeeze(imag(GAF_loew(kk,jj,:))),linestyles{ii},'LineWidth',2)
        end
    end
    % Check eigs.
    % if (max(real(eig(sys_loew))) > 0)
    %     'Unstable system'
    % end
end
legend('Ref', [int2str(Nmax(1)) ' states'],...
              [int2str(Nmax(2)) ' states'],...
              [int2str(Nmax(3)) ' states'])


% eof


