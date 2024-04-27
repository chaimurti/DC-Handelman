%This function takes two points in the x,y plane and returns the
%a vector of the coefficients of [y x 1] of the line connecting the two points
%
%Details
%
%Inputs:
%A,B := 1 x 2 vectors
%
%Outputs:
%out:= 1 x 3 vector
%
%%

function out = handelman_linemaker_2d(A,B);

out = [A(1)*(B(2)-A(2))-A(2)*(B(1)-A(1)), -(B(2)-A(2)), B(1)-A(1) ];
if out(1)~=0
    out = (1/abs(out(1)))*out;
else
    out = out;
end

 


