% Function Sears. 
%  Compute the analytical expression for Sears's function about the 
%  leading edge in terms of modified Bessel functions of the second kind.
%
%  Inputs: an array of reduced frequencies, k_n, as real numbers.
%  Outputs: complex values of S_0(i k_n).
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S0=sears(k)
  S0=zeros(size(k));
  for n=1:length(k)
     if k(n)==0
         S0(n)=1;
     elseif k(n) > 0
         S0(n)= exp(-1i*k(n)) ...
              ./(1i*k(n).*(besselk(0,1i*k(n))+besselk(1,1i*k(n))) );
     else
         error('Input must be non-negative.')
     end
  end
end