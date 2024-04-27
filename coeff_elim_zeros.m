% This function takes as its inputs coefficient matrix ('coeffs'),
% a matrix of exponents ('exps'), and a 'facet_rule'. It returns a matrix of coefficients
%where all columns corresponding to elements of 'exps' where the only nonzero elements correspond to 
%'facet_rule' are removed.
%
%Details:
%Inputs: 
%Coeffs:= m x n matrix
%exps:= m x d matrix
%facet_rule:= 1 x d vector
%
%Outputs:
%out_coeffs:= m x d' (d' < d) matrix
%
%%

function [out_coeffs] = coeff_elim_zeros(Coeffs,exps,facet_rule)

Kk = find(facet_rule);
S = size(Coeffs);
S = S(1);
out_coeffs = zeros(S,0);
J=1;
T = size(exps);
length = T(1);
temp = zeros(S,1);
for i=1:length
    tempexp = exps(i,:);
    test=sum(tempexp(Kk));
    if test~=0
        out_coeffs(:,i) = temp;
    else
        out_coeffs(:,i) = Coeffs(:,i);;
    end
end