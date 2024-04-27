% This function takes as its inputs coefficient matrix ('coeffs'),
% a matrix of exponents ('exps'), and a 'facet_rule'. It returns a matrix of coefficients
%where all columns corresponding to elements of 'exps' where the only nonzero elements correspond to 
%'facet_rule' are removed. Also returned is the modified vector of exponents
% and the total number of coefficients after the modification. This
% performs the work of the map R_i(b_i,d)
%
%Details:
%Inputs: 
%Coeffs:= m x n matrix
%exps:= m x d matrix
%facet_rule:= 1 x d vector
%
%Outputs:
%out_coeffs:= m x d' (d' < d) matrix
%out_exps:= m x d matrix
%Coeff_total:= m
%
%%
function [out_coeffs,out_exps,out_coeff_total] = coeff_elim(Coeffs,exps,facet_rule)

Kk = find(facet_rule);
S = size(Coeffs);
S = S(1);
out_coeffs = zeros(S,0);
J=1;
T = size(exps);
length = T(1);
T = T(2);
out_exps = zeros(0,T);
for i=1:length
    tempexp = exps(i,:);
    test=sum(tempexp(Kk));
    if test~=0
        out_coeffs(:,J) = Coeffs(:,i);
        out_exps(J,:) = exps(i,:);
        J=J+1;
    end
end

TT = size(out_coeffs);
out_coeff_total = TT(2);