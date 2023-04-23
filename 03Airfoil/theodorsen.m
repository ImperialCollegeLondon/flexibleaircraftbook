% Function theodorsen. 
%  Compute the analytical expression for Theodorsen's lift deficiency
%  function in terms of modified Bessel functions of the second kind.
%
%  Inputs: an array of reduced frequencies, k_n, as real numbers.
%  Outputs: complex values of C(i k_n).
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: April 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=theodorsen(k)
  C=zeros(size(k));
  for n=1:length(k)
     if k(n)==0
         C(n)=1;
     elseif k(n) > 0
         C(n)=besselk(1,1i*k(n))./(besselk(0,1i*k(n))+besselk(1,1i*k(n)));
     else
         error('Input must be non-negative.')
     end
  end
end
