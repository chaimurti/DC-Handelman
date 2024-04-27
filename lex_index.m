%This function returns the reverse lexicographical index of the monomial
%with exponents 'alfa' of a homogeneous polynomial of degree 'd' in 'n' variables, 
%
%Details
%
%Inputs:
%n:= positive scalar
%d:= positive scalar
%alfa:= 1 x n vector
%
%Outputs:
%I:= positive scalar
%
%%

function I=lex_index(alfa,n,d)
I=0;
for j=1:n-1
    J=0;
    s=0;
    for k=1:j-1
        s=s+alfa(k);
    end
    for i=1:alfa(j)
        J=J+factorial(n-j+d+1-s-i-1)/(factorial(d+1-s-i)*factorial(n-j-1));
    end
    I=I+J;
end
I=I+1;