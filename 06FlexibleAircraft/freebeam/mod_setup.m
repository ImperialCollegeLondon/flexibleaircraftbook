%% Compute problem matrices.
% Assemble stiffness matrix.
Kfull=sparse(2*(N+1),2*(N+1));
for i=1:N
  ldofs=2*i-1:2*(i+1); % Local degrees of freedom
  Kfull(ldofs,ldofs)= Kfull(ldofs,ldofs) + Kmat(EI,L/N);
end
Kss=Kfull(adofs,adofs); % Apply boundary conditions

K0=zeros(2*N+3);
K0(4:end,4:end)=Kss;


% Compute constant terms in mass matrix
m=rho_s*L;
rcg=[L/2; 0];
rcgcross=[-rcg(2);rcg(1)];
IB=(m/3)*L^2;
Mrx_full=zeros(2,2*(N+1));
Mpx_full=zeros(1,2*(N+1));
Mxx_full=sparse(2*(N+1),2*(N+1));
for i=1:N
  ldofs=2*i-1:2*(i+1); % Local degrees of freedom
  Mrx_full(2,ldofs)    = Mrx_full(2,ldofs)     + Mrx_e(rho_s,L/N);
  Mpx_full(1,ldofs)    = Mpx_full(1,ldofs)     + Mpx_e(rho_s,rBP(i:i+1),L/N);
  Mxx_full(ldofs,ldofs)= Mxx_full(ldofs,ldofs) + Mxx_e(rho_s,L/N);
end
Mrx=Mrx_full(1:2,adofs);
Mpx=Mpx_full(1,adofs);
Mxx=Mxx_full(adofs,adofs);
clear Mrx_full Mpx_full

M0=[[m*eye(2)      -m*rcgcross Mrx]; ...
    [-m*rcgcross'     IB       Mpx]; ...
    [Mrx'             Mpx'     Mxx]];
invM0=inv(M0);
invMrigid=inv(M0(1:3,1:3));

% eof