function out=post_snaps(r,v,w,dt,Nframe)

    % Compute instantaneous orientation of reference frame.
    theta=-cumtrapz(v(:,3))*dt;
    
    % Compute instantaneous inertia velocities of reference frame.
    vx=v(:,1).*cos(theta)-v(:,2).*sin(theta);
    vz=v(:,1).*sin(theta)+v(:,2).*cos(theta);

    % Compute instantaneous inertial coordinates of reference frame.
    ux0(:)=cumtrapz(vx)*dt;
    uz0(:)=cumtrapz(vz)*dt;

    Nt=length(vx);
    
    % Compute snapshots and displace them by DeltaX.
    DeltaX=r(end)/8;
    
    figure
    k=0;
    for i=1:Nframe:Nt
        k=k+1;
        ux(1)=ux0(i)+k*DeltaX;
        uz(1)=uz0(i);
        for j=1:length(r)-1
            ux(j+1)=ux(1)+r(j+1)*cos(theta(i))-w(i,j)*sin(theta(i));
            uz(j+1)=uz(1)+r(j+1)*sin(theta(i))+w(i,j)*cos(theta(i));
        end
       plot(ux,uz,'k'), hold on
    end
    %axis equal
    %axis tight
    

    % Plot the inertial velocities at the origin.
    if 0
        t=[0:dt:(Nt-1)*dt];   
        figure
        subplot(2,1,2)
            plot(t,theta*180/pi)
            xlabel('time')
            ylabel('\theta (deg)')

        subplot(2,1,1)
            plot(t,vx,'k--','LineWidth',2), hold on
            plot(t,vz,'k','LineWidth',2)
            legend('Vx','Vz')
    end
    out.theta=theta;
    out.vx=vx;
    out.vz=vz;
end