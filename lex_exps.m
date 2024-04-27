%This file returns a matrix of exponents of the monomials of a
%non-homogeneous polynomial in 'var' indetermintes and of degree 'deg'.
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
function out = lex_exps(vars,deg)

Vec = zeros(1,vars);
for j=1:deg
    dummy = homopoly(vars,j);
    Vec = [Vec;dummy];
end
d = size(Vec);
d = d(1);
for j=1:d
    k = lex_index_nh(Vec(j,:));
    out(k,:) = Vec(j,:);
end