%This function returns the total number of Handelman coefficients of a polynomial of degree 'deg'
%on a polytope with 'n_polytope' facets.
%Details
%
%Inputs:
%n_polytope:= positive scalar > 2
%deg:=  positive scalar
%
%Outputs:
%out:= positive scalar

function out=Coeff_total(n_polytope,deg)


out = 0;

for i=0:deg
    dummy = (factorial(i+n_polytope-1))/(factorial(i)*factorial(n_polytope-1));
    out = out+dummy;
end
