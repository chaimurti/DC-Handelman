%This function eliminates redundant constraints: that is, if 
%two rows of a constraint matrix are equivalent, one is removed.
%
%Details
%
%Inputs:
%A:= m x n matrix
%
%Outputs
%Aout:= m' x n matrix, m' < m
%
%%
function  Aout = constraint_elim(A)

size_A = size(A);
rows_A = size_A(1);
Aout = [];
k = 1;
zerosmat = zeros(1,size_A(2));

for i=1:size_A(1)-1
    if sum(find(A(i,:)))>0
        for j=i+1:size_A(1)
            if sum(find(A(i,:)-A(j,:)))==0
                A(j,:) = zerosmat;
            end
            if sum(find(A(i,:)+A(j,:)))==0
                A(j,:) = zerosmat;
            end
        end
    end
end


% 
for i=1:size_A(1)
    if sum(find(A(i,:)))
        Aout(k,:) = A(i,:);
        k = k+1;
    end
end
