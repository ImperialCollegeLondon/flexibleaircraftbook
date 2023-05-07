% Function theodorsen_rfa. 
%  Compute a rational-function approximation to Theodorsen's lift 
%  deficiency function of order Na.
%
%  Inputs: Order of the approximation (Na).
%  Outputs: State-space object for the SISO LTI approximing C(ik).
%
% Dependencies:
%    theodorsen.m: Analytical expression for Theodorsen's lift deficiency
%                  function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: May 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sysC=theodorsen_rfa(Na)

    % We obtain the RFA via sampling of the FRF. Define reduced frequencies
    % at which the function will bé sampled.
    % Note: C(ik) is a rather smooth function and results are pretty 
    % robust to this choice.
    deltak=0.01;
    k=0:deltak:5;
    k=[k 100];

    % Add extra weight to the samples for k<1
    Wk=ones(size(k));
    Wk(1:1/deltak)=100;

    % Define the data-driven FRF, removing the feedthrough contribution of
    % C(ik) to speed up convergence in the fitting process, and obtain
    % the best fit state-space model using a log-Chebyshev approximation.
    systh=frd(theodorsen(k)-0.5,k);
    sysC=fitmagfrd(systh,Na,1,Wk)+0.5;

end

% eof