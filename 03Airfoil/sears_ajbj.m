% Function sears_ajbj. 
%  Compute the coefficients a_j and b_j that result when the RFA to 
%  Sears's function is written in the form 1-sum_j((aj*x)(bj+x).
%
%  Inputs: Order of the approximation (Na).
%  Outputs: a_j and b_j coefficients approximating S_0(ik).
%
% Dependencies:
%    sears_rfa.m: RFA to Sears's function.
%
% Written by: Rafael Palacios (r.palacios@imperial.ac.uk)
% Latest update: May 2023. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,b]=sears_ajbj(Na)

    % First, obtain the RFA to Sears of order Na
    sysS0=sears_rfa(Na);

    % We rewrite this expression as 1-sum_j((aj*x)(bj+x). The new 
    % coefficients result from a linear system from s_0(ik) written in 
    % terms of poles and zeros by equating both expressions. 
    % This is done next usign symbolic algebra.
    
    % Next, the b_j coefficients are the poles of the SS obtained above.
    [z,p,k]=zpkdata(sysS0);
    b=-p{1};

    % Polynomials are written as a function of x.
    syms x

    % Define the unknowns as symbolic variables.
    aj=sym('a',[1 Na]); 
    bj=sym('b',[1 Na]);  % We know this one, but we keep symbolic for now.

    % Two matrices to store the polynomials in the problem.
    C1=1;  % C1 stores the RFA in terms of aj and bj.
    C2=1;  % C2 stores the denominator of the zpk form.
    for j=1:Na
        C1=C1-aj(j)*x/(bj(j)+x);
        C2=C2*(bj(j)+x);
    end

    % Compute the coefficients of the polynomial in x, with known bj. 
    CC=subs(coeffs(collect(C1*C2,x),x),bj,b');

    % Each coefficient  has to be equal to the numerator of the 
    % zpk form of S0(ik), which gives linear equations in aj 
    eqs=k*poly(z{1}) - CC(end:-1:1);
    
    % Solve the equation in the a_j coefficients.
    [AA,bb]=equationsToMatrix(eqs(1:Na),aj);
    a=double(AA\bb);
end
% eof