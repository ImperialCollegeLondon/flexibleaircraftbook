%SHAPEDER   Derivative of the shape function. 
%   SHAPEDER(x,L,n) computes the derivative of the shape function at point 
%   x within a bending element of length L and 4 degrees of freedom . 
%   The degrees of freedom are 
%     n=1: displacement at first node; 
%     n=2: rotation at first node;
%     n=3: displacement at second node;
%     n=4: rotation at second node.   
function Nder=shapeder(x,L,n)

if any(x<0) || any(x>L)
  disp('x out of bounds in SHAPE routine')
  return
else
  if n==1
    Nder=   - 6*x/L^2 + 6*x.^2/L^3;
  elseif n==2
    Nder= 1 - 4*x/L   + 3*x.^2/L^2;
  elseif n==3
    Nder=     6*x/L^2 - 6*x.^2/L^3;
  elseif n==4
    Nder=   - 2*x/L   + 3*x.^2/L^2;
  end
end

% eof