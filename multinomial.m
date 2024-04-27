%This function provides the expansion of a monomial with coefficients 
%in 'vec' raised to the power 'deg'
%
%Details
%
%Inputs:
%vec:= 1 x n vector
%deg:= positive scalar
%
%Outputs:
%out:= 1 x m vector, m > n if deg >= 2
%
%%
function out = multinomial(vec,deg)

if deg==0
    out = 1;
else
    vars = length(vec);

temp_exps = homopoly(vars,deg);
temp_exps2 = zeros(size(temp_exps));
J = length(temp_exps);
for i=1:J
    temp_exps2(i,:) = temp_exps(J+1-i,:);
end
temp_exps = temp_exps2;
temp_out = zeros(1,J);
temp_prod_fac = 1;
temp_prod_coeff=1;
for i=1:J
    for j=1:length(vec)
        temp_prod_fac = temp_prod_fac*factorial(temp_exps(i,j));
        temp_prod_coeff = temp_prod_coeff*vec(j)^(temp_exps(i,j));
    end
    k = lex_index_nh(temp_exps(i,2:vars));
    temp_out(k) = (factorial(deg)/temp_prod_fac)*temp_prod_coeff;
    temp_prod_fac = 1;
    temp_prod_coeff=1;
end
out = temp_out;
end
    
