%This function gives the lex index of a monomial with exponent 'alfa' in some non-homogeneous polynomial.
%The degree and number of variables of the polynomial is determined internally.
%
%Details
%
%Inputs:
%alfa:= 1 x n vector
%
%Outputs:
%lexnh:= positive scalar
%
%%
function lexnh = lex_index_nh(alfa)

deg = sum(alfa);
n = length(alfa);

if deg==0
    lexnh = lex_index(alfa,n,0);
else
    lexnh = Coeff_total(n,deg-1) + lex_index(alfa,n,deg);
end
