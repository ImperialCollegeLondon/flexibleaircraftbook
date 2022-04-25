function Mpx=Mpx_e(rho,rBP,L)

for j=1:4
    Mfun = @(x)-rho*(rBP(1)+(x/L)*(rBP(2)-rBP(1))).*shape(x,L,j);
    Mpx(j)=quadl(Mfun,0,L);
end

% eof