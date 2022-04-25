% Compute linear normal modes from solution in quasi-coordinates. 
[V1,D1]=eig(full(K0),full(M0));
omega=sqrt(diag(D1));

%Mass normalization (pick structural modes with positive tip deflection)
for j=1:2*N+3
    if j>3
        factor=sign(V1(end-1,j));
    else
        factor=1;
    end
    norm=transpose(V1(:,j))*full(M0)*V1(:,j)*factor;
    V1(:,j)=V1(:,j)/norm;
end


% Define block matrices.
PhiRR=V1(1:3,1:3);
PhiRS=V1(1:3,4:end);
PhiSS=V1(4:end,4:end);


% Plot the nodal vertical displacements for the first four modes.
figure
subplot(2,2,1)
 plot(0:1/N:1,[0; PhiSS(1:2:end,1)],'k','LineWidth',2)
ylabel('\Phi_{SS,1}')
subplot(2,2,2)
 plot(0:1/N:1,[0; PhiSS(1:2:end,2)],'k','LineWidth',2)
ylabel('\Phi_{SS,2}')
subplot(2,2,3)
 plot(0:1/N:1,[0; PhiSS(1:2:end,3)],'k','LineWidth',2)
xlabel('x/L')
ylabel('\Phi_{SS,3}')
subplot(2,2,4)
 plot(0:1/N:1,[0; PhiSS(1:2:end,4)],'k','LineWidth',2)
xlabel('x/L')
ylabel('\Phi_{SS,4}')

% Show integral relations between PhiSS and PhiRS.
disp(['Check relation between PhiSS and PhiRS = ' ...
num2str(max(abs(invMrigid*[Mrx;Mpx]*PhiSS+PhiRS),[],'all'))])


% Transform mode shapes into absolute displacements.
for j=2:2*N+3
  V1B(1:2:2*(N+1),j)=[0; V1(4:2:end,j)] + V1(2,j)*ones(N+1,1) + ...
                    -(V1(3,j)*[0:L/N:L]');
end
figure
subplot(3,2,1)
 plot(0:1/N:1,V1B(1:2:end,2),'k','LineWidth',2)
 ylabel('Mode 2'), axis([0 1 -1 1])
subplot(3,2,2)
 plot(0:1/N:1,V1B(1:2:end,3),'k','LineWidth',2)
 ylabel('Mode 3'), axis([0 1 -1 1])
subplot(3,2,3)
 plot(0:1/N:1,V1B(1:2:end,4),'k','LineWidth',2)
 ylabel('Mode 4'), axis([0 1 -1 1])
subplot(3,2,4)
 plot(0:1/N:1,V1B(1:2:end,5),'k','LineWidth',2)
 ylabel('Mode 5'), axis([0 1 -1 1])
subplot(3,2,5)
 plot(0:1/N:1,V1B(1:2:end,6),'k','LineWidth',2)
 ylabel('Mode 6'), axis([0 1 -1 1]), xlabel('x/L') 
subplot(3,2,6)
 plot(0:1/N:1,V1B(1:2:end,7),'k','LineWidth',2)
 ylabel('Mode 7'), axis([0 1 -1 1]), xlabel('x/L')

%% Check: Compute modes from finite-element matrices directly.
if 0
    [V2,D2]=eig(full(Kfull),full(Mxx_full));
    % Mass normalization
     for j=1:2*(N+1)
         norm=transpose(V2(:,j))*full(Mxx_full)*V2(:,j);
         V2(:,j)=-V2(:,j)/norm;
     end

    figure
    subplot(3,2,1)
     plot(0:1/N:1,V2(1:2:end,1)),  axis([0 1 -1 1])
    subplot(3,2,2)
     plot(0:1/N:1,V2(1:2:end,2)),  axis([0 1 -1 1])
    subplot(3,2,3)
     plot(0:1/N:1,V2(1:2:end,3)),   axis([0 1 -1 1])
    subplot(3,2,4)
     plot(0:1/N:1,V2(1:2:end,4)),   axis([0 1 -1 1])
    subplot(3,2,5)
     plot(0:1/N:1,V2(1:2:end,5)),   axis([0 1 -1 1])
    subplot(3,2,6)
     plot(0:1/N:1,V2(1:2:end,6)),   axis([0 1 -1 1])
end
 

% eof