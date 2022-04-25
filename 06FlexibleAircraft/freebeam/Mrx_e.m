function Mrx=Mrx_e(rho,L)

for j=1:4
    Mfun = @(x)rho*shape(x,L,j);
    Mrx(j)=quadl(Mfun,0,L);
end

% eof