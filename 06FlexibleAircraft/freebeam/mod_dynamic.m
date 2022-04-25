% Initial conditions:
x0 = zeros(4*N+3,1);

% Rotation with constant angular velocity at the centre of mass.
% x0(3,1)=0.1*pi/2;  
% x0(2,1)= x0(3,1)*L/2;

% Constant normal force at the free end.
F=zeros(2*N+3,1);
F(3,1)=-Fmax*L;
F(3+2*N-1,1)=Fmax;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve rigid problem.
% For pure moment, solution is a omega=at^2 for t<Tmax, with a=FmaxL/(2TmaxIB_cg)
 [t,x_rigid]=ode45(@frigid,[0:dt:tfinal],x0(1:3,1));
 
% out_rigid=post_snaps(rBP,x_rigid,zeros(length(t),N),dt,10);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modal solution
xmodal0=zeros(3+2*Nflexmodes,1);
Fmodal=transpose(V1(:,1:3+Nflexmodes))*F;
[t,xmodal]=ode45(@flinear,[0:dt:tfinal],xmodal0);
x1=[xmodal(:,1:3) xmodal(:,Nflexmodes+4:end)]*transpose(V1(:,1:3+Nflexmodes));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Plot modal amplitudes.
figure
subplot(3,1,1)
plot(t,xmodal(:,3),'k','LineWidth',2)
ylabel('q_3','FontSize',16)
legend(['EI=' int2str(EI) ' Nm^2'])
subplot(3,1,2)
plot(t,xmodal(:,Nflexmodes+5),'k','LineWidth',2)
ylabel('q_5','FontSize',16)
subplot(3,1,3)
plot(t,xmodal(:,Nflexmodes+7),'k--','LineWidth',2)
ylabel('q_7','FontSize',16)
xlabel('t (s)','FontSize',16)



% Create beam deformation plot.
out1=post_snaps(rBP,x1(:,1:3),x1(:,4:2:end),dt,20);

% Compare rigid and modal solution at the free end.
figure
subplot(1,2,1)
    plot(t,x_rigid(:,2),'k:','LineWidth',2)
    hold on
    plot(t,x1(:,2),'k--','LineWidth',2)
    ylabel('v_z (m/s)')
    legend('Rigid',['EI=' int2str(EI) ' Nm^2'])
    grid on

subplot(1,2,2)
    plot(t,x_rigid(:,3),'k:','LineWidth',2)
    hold on
    plot(t,x1(:,3),'k--','LineWidth',2)
    ylabel('\omega_y (rad/s)')
    grid on

    
figure    
    plot(t,x1(:,end-1)/L,'k--','LineWidth',2)
    ylabel('w(L)/L')
    xlabel('t (s)')
    grid on
    




figure
subplot(2,1,1)
    plot(t,x1(:,1:3))
    grid on
    legend('Vx','Vz','Wy')

subplot(2,1,2)
    plot(t,x1(:,4:2:end)/L)
    ylabel('w/L')
    xlabel('time')
    grid on

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem using stiff ODE solver.

[t,x2]=ode15s(@fnln,[0:dt:tfinal],x0);

figure
subplot(2,1,1)
    plot(t,x2(:,1:3))
    grid on
    legend('Vx','Vz','Wy')

subplot(2,1,2)
    plot(t,x2(:,2*N+4:2:end)/L)
    ylabel('w/L')
    xlabel('time')
    grid on

figure
    for i=1:length(t)
     KinEnergy(i)=1/2* x2(i,1:2*N+3)*M0*x2(i,1:2*N+3)';
     PotEnergy(i)=1/2* x2(i,2*N+4:end)*Kss*x2(i,2*N+4:end)';
    end
    plot(t,KinEnergy,'k--'), hold on
    plot(t,PotEnergy,'k')
    plot(t,KinEnergy+PotEnergy)
    grid on

out2=post_snaps(rBP,x2(:,1:3),x2(:,2*N+4:2:end-1),dt,10);

%post_video(VideoWriter('flyingbar.avi'),...
%           rBP,out2,x2(:,2*N+4:2:end-1),dt,size(t)/100);

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function dxdt=flinear(t,x)
    global Nflexmodes Tmax omega Fmodal
    A=zeros(2*Nflexmodes+3);
    for j=1:Nflexmodes
        A(3+j,3+Nflexmodes+j)=-omega(j+3)^2;
        A(3+Nflexmodes+j,3+j)=1;
    end
    Force=[Fmodal(1:Nflexmodes+3,1); zeros(Nflexmodes,1)];
    if t<Tmax
      Force=Force*t/Tmax;
    end
    dxdt = A * x + Force;
end


% Nonlinear rigid solution.    
function dxdt = frigid(t,x)
  global invM0 Kss m IB rcgcross Mrx Mpx Mxx F N
  global Tmax
  global invMrigid
    vB=x(1:2);
    oB=x(3);
    
    oBcross=[[0 oB]; [-oB 0]];

    C=[[m*oBcross       -m*oBcross*rcgcross];
       [zeros(1,2)       m*rcgcross(2)*vB(1)]];
   
    Force=F(1:3,1);
    if t<Tmax
      Force=Force*t/Tmax;
    end
   dxdt = invMrigid*(Force - C*x);
end


% Nonlinear solution.    
function dxdt = fnln(t,x)
  global invM0 Kss m IB rcgcross Mrx Mpx Mxx F N
  global Tmax
    Nxi=(length(x)-3)/2;
    vB=x(1:2);
    oB=x(3);
    oBcross=[[0 oB]; [-oB 0]];

    C=zeros(Nxi+3);
    C=[[m*oBcross       -m*oBcross*rcgcross    2*oBcross*Mrx];
       [zeros(1,2)       m*rcgcross(2)*vB(1)   2*oB*Mpx];
       [Mrx'*oBcross     -oB*Mpx'             zeros(Nxi)]];
       
    A=[[-C [zeros(3,Nxi); -Kss]];
       [zeros(Nxi,3) eye(Nxi) zeros(Nxi)]];
   
    Force=[F; zeros(Nxi,1)];
    if t<Tmax
      Force=Force*t/Tmax;
    end
   dxdt =[[invM0    zeros(Nxi+3,Nxi)]; ...
          [zeros(Nxi,Nxi+3) eye(Nxi)  ]] ...
         * (A * x + Force);
end
% eof