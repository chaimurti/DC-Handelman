%This function takes two polynomials of 'vars' variables, with 
%coefficients in 'vecs' and with degrees in 'degs', and returns a 
%vector of coefficients of their product
%
%Details
%
%Inputs:
%vecs:= 1 x k cell 
%vars:= positive scalar
%degs:= 1 x k vector
%
%Outputs
%out:= 1 x m vector
%
%%
function out = polyprod2(vecs,vars,degs)

deg_total = sum(degs);
temp_coeff_total = Coeff_total(vars,deg_total);
temp_out = zeros(1,temp_coeff_total);
K0 = size(vecs);
num_vecs = K0(1,1);
outlength=length(vecs{1});
test2 = vecs{1};
temp_deg = degs(1);
out_exps = lex_exps(vars,temp_deg);

for i=2:num_vecs
    temp_vec = vecs{i};
    temp_vec_exps = lex_exps(vars,degs(i));
        for k=1:outlength;
            test = test2(k)*temp_vec(1);
            testexp =temp_vec_exps(1,:)+out_exps(k,:);
            l = lex_index_nh(testexp);
            temp_out(l)= test;           
        end 
        for j=2:length(temp_vec)         
                for k=1:outlength;
                    test = test2(k)*temp_vec(j);
                    testexp =temp_vec_exps(j)+out_exps(k);
                    l = lex_index_nh(testexp);
                    temp_out(l)= temp_out(l) + test;           
                end   

        end
    temp_deg = temp_deg+degs(i);
    out_exps = lex_exps(vars,temp_deg);
    outlength = Coeff_total(vars,temp_deg);
    test2 = temp_out(1:outlength);
end
out = temp_out;

                  