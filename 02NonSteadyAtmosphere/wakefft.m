% Function discrete_fft
%
%  Compute the frequency content (in wavenumber) of 1-cos and vertical 
%  wake vortex profiles.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: January 2018
%%%%%%%%%%%%%%%%%%%
clear all, close all

Fs = 1000;          % Sampling spatial frequency
t = 0:1/Fs:50;      % Vector of coordinates.
L=length(t);

% Loop through core values.
rc=[0.03 0.05 0.1]
leg={'r_c=0.03b',...
     'r_c=0.05b',...
     'r_c=0.1b'};
lines={'k-','k--','k:'};


for j=1:length(rc)
% Wake -- coordinates normalized by b.
        t0=max(t)/2;
    for k=1:L
        r1=-(t(k)-t0);
        r2= t(k)-(t0+pi/4);
        X(k)=r1/(r1^2+rc(j)^2)+r2/(r2^2+rc(j)^2);
    end
    %

    figure(1)
    plot(t,X,lines{j}), hold on
    xlabel('x/(b)')
    ylabel('w_g(t)')

    % To use the fft function to convert the signal to the frequency domain,
    % first identify a new input length that is the next power of 2 from the 
    % original signal length. This will pad the signal X with trailing zeros 
    % in order to improve the performance of fft.
    n = 2^nextpow2(L);

    % Convert the given signal to the frequency domain.
    Y = fft(X,n);

    % Define the frequency domain and plot the single-sided wavenumber spectrum.
    f = Fs*(0:(n/2))/n;
    P = abs(Y/n);
    P = P/max(P); % normalize

    figure(2)
    plot(f,P(1:n/2+1),lines{j},'LineWidth',2), hold on
end
xlabel('\kappa b/(2\pi)','FontSize',14)
ylabel('FFT(w_g)','FontSize',14)
axis([0 10 0 1])
legend(leg)