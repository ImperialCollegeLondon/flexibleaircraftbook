%SHAPE   Shape function. 
%   SHAPE(x,L,n) computes the shape function at point x within a bending
%   element of length L and 4 degrees of freedom. The degrees of freedom
%   are 
%     n=1: displacement at first node; 
%     n=2: rotation at first node;
%     n=3: displacement at second node;
%     n=4: rotation at second node.   
function Nfun=shape(x,L,n)

if any(x<0) || any(x>L)
  disp('x out of bounds in SHAPE routine')
  return
else
  if n==1
    Nfun= 1 - 3*x.^2/L^2 + 2*x.^3/L^3;
  elseif n==2
    Nfun= x - 2*x.^2/L   +   x.^3/L^2;
  elseif n==3
    Nfun=     3*x.^2/L^2 - 2*x.^3/L^3;
  elseif n==4
    Nfun=   -   x.^2/L   +   x.^3/L^2;
  end
end

