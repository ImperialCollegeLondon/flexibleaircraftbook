%KMAT Material stiffness matrix. 
%   Kmat (EI,L) computes the material stiffness matrix for a beam
%   bending element of length L, bending stiffness EI. The beam has 4 
%   degrees of freedom, given by:
%     i,j=1: displacement at first node; 
%     i,j=2: rotation at first node;
%     i,j=3: displacement at second node;
%     i,j=4: rotation at second node.
function StiffMat=Kmat(EI,L)

for i=1:4
  for j=1:4
    KFun = @(x)EI*shapeder2(x,L,i).*shapeder2(x,L,j);
    StiffMat(i,j)=quadl(KFun,0,L);
  end
end

% eof