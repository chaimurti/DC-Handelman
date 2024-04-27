clear model;
clear params;
A=[-A_ineq ; Aeq];
model.A = sparse(A);
model.obj = [zeros(1,length(cmat))];
model.modelsense = 'Min';
model.rhs = [-b_ineq', bmat'];

for (i=1:length(cmat))
    ineq(i)='<';
end
for (i=1:length(bmat))
    eq(i)='=';
end

model.sense=[ineq,eq];
params.OptimalityTol=1e-9;
params.FeasibilityTol=1e-9;
params.ScaleFlag=1;
% params.BarCorrectors=10;
params.Method=2;
result = gurobi(model,params);
result
max(result.x)
min(result.x)


% clear params;
% params.method = 2;
% params.timelimit = 100;
% 
% result2 = gurobi(model, params);
% result2
% % disp(result.objval);
% % disp(result.x)
% max(result2.x)
% min(result2.x)
