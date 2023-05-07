% Function sears_rfa. 
%  Compute a rational-function approximation to Sears's function of order Na.
%
%  Inputs: Order of the approximation (Na).
%  Outputs: State-space object for the SISO LTI approximing S_0(ik).
%
% Dependencies:
%    sears.m: Analytical expression for Sears's function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: May 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sysS0=sears_rfa(Na)

    % We obtain the RFA via sampling of the FRF. Define reduced frequencies
    % at which the function will bé sampled. A large value of k_max is needed to 
    % fit the asymptotic value through fitmagfrd (which does
    % enforce that S goes to zero automatically).
    deltak=0.005; 
    k=[[0:deltak:5-deltak] [5:deltak*20:50]];

    % Extra weights to points for k<1, where potential-flow captures the 
    % physics.
    Wk=ones(size(k));
    Wk(1:1/deltak)=100;
    
    sysse=frd(sears(k),k);
    sysS0=fitmagfrd(sysse,Na,[],Wk);

end

% eof