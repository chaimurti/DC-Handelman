%This file returns a matrix of exponents of the monomials of a homogeneous
%polynomial in 'var' indetermintes and of degree 'deg'.
%
%Details
%
%Inputs:
%var:= positive scalar
%deg:= positive scalar
%
%Outputs:
%exponents:= m x var matrix
%
%%


function exponents = homopoly(var,deg)
exponents = [];
if (var == 1)
    exponents = deg;
    return
elseif (deg == 0)
    exponents = zeros(1,var);
    return
end
for i=deg:-1:0
    dummy = homopoly(var-1,deg-i);
    exponents = [exponents; i*ones(size(dummy,1),1) dummy];
end
end