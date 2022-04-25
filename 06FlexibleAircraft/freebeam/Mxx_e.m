function Mxx=Mmat(rho,L)

for i=1:4
  for j=1:4
    Mfun = @(x)rho*shape(x,L,i).*shape(x,L,j);
    Mxx(i,j)=quadl(Mfun,0,L);
  end
end

% eof