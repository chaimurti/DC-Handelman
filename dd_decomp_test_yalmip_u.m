
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Local DC Decomposition of nonconvex polynomials using Handelman's Theorem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The main objective of this script is to accept a polynomial P(x), and
%returna two convex polynomials gamma(x) and nu(x) such that P = gamma-nu
%and both gamma and nu are convex. We do so using LP/SDP relaxations after
%applying Handelman's Theorem.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initializing state space, P(x), gamma(x), nu(x), Gamma_1, and Gamma_2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

display('Initializing state space, P(x), gamma(x), nu(x), Gamma_1, and Gamma_2')
relax_type = 'SDP';    
n = 30;                                                  %Dimension of state space





P.deg = 4;                                              %Degree of polynomial to be decomposed
offset_deg = 0;                                         %Degree offset
c1 = -5;
% c2 = 10;

gamma.deg = P.deg+offset_deg;                           %Degree of gamma(x)
nu.deg = P.deg+offset_deg;                              %Degree of nu(x)
P.monoms = lex_exps(n,P.deg);                           %Monomials of P(x)
% P.coeffs = zeros(1,Coeff_total(n,P.deg));   
P.coeffs = 7*(rand(1,Coeff_total(n,P.deg))-0.5); %Set up a vector of coefficients for P(x)
P.coeff_total = Coeff_total(n,P.deg);
% temp_lex = lex_index_nh([4 0]);
% P.coeffs(temp_lex) = 1;
% temp_lex = lex_index_nh([1 3]);
% P.coeffs(temp_lex) = -1;
% temp_lex = lex_index_nh([2 1]);
% P.coeffs(temp_lex) = -1;
% temp_lex = lex_index_nh([1 1]);
% P.coeffs(temp_lex) = c1;
% temp_lex = lex_index_nh([2 0]);
% P.coeffs(temp_lex) = -1;
% 
% temp_lex = lex_index_nh([1 0]);
% P.coeffs(temp_lex) = -2;
% temp_lex = lex_index_nh([0 1]);
% P.coeffs(temp_lex) = 3;
% temp_lex = lex_index_nh([0 0]);
% P.coeffs(temp_lex) = 1;

%Degree offset

gamma.deg = P.deg+offset_deg;                           %Degree of gamma(x)
nu.deg = P.deg+offset_deg;                              %Degree of nu(x)



gamma.monoms = lex_exps(n,gamma.deg);                   %Monomials of gamma(x)
gamma.coeff_total = Coeff_total(n,gamma.deg);            %Total number of coefficients in gamma(x)
gamma.coeffs = ones(1,Coeff_total(n,gamma.deg));        %Set up a vector of coeffs for gamma(x); these are decision variables
gamma.constraint_1 = diag(gamma.coeffs);                %Set up a constraint matrix for gamma(x)

nu.monoms = lex_exps(n,nu.deg);                         %Monomials of nu(x)
nu.coeff_total = Coeff_total(n,nu.deg);                 %Total number of coefficiens in nu(x) 
nu.coeffs = ones(1,Coeff_total(n,nu.deg));              %Set up a vector of coeffs for nu(x); these are decision variables
nu.constraint_1 = diag(nu.coeffs);                      %Set up a constraint matrix for gamma(x)


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating Hessian matrices for gamma(x) and nu(x).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating Hessian matrices for gamma(x) and nu(x)...')

gammma.grad_monoms = zeros(2,Coeff_total(n,gamma.deg-1));       %Set up monomials for gradient of gamma(x).
nu.grad_monoms = zeros(2,Coeff_total(n,nu.deg-1))';             %Set up  monomials for gradient of nu(x).

gammma.hess_monoms = zeros(2,Coeff_total(n,gamma.deg-2));       %Set up monomials for gradient of gamma(x).
nu.hess_monoms = zeros(2,Coeff_total(n,nu.deg-2))';             %Set up monomials for gradient of nu(x).

%Set up gradient of gamma(x)
for i=1:n
    temp = gamma.monoms(:,i) - ones(Coeff_total(n,gamma.deg),1);
    temp_coeffs = gamma.monoms(:,i).*gamma.coeffs(:);
    temp_constr = gamma.constraint_1*diag(gamma.monoms(:,i));
    gamma.grad_monoms(:,i) = temp(temp>=0);
    gamma.grad_coeffs(:,i) = temp_coeffs(temp_coeffs>0);
    gamma.grad_constraint(:,:,i) = temp_constr(any(temp_constr,2),:);    
end


%Set up gradient of nu(x)
for i=1:n
    temp = nu.monoms(:,i) - ones(Coeff_total(n,nu.deg),1);
    temp_coeffs = nu.monoms(:,i).*nu.coeffs(:);
    temp_constr = nu.constraint_1*diag(nu.monoms(:,i));
    nu.grad_monoms(:,i) = temp(temp>=0);
    nu.grad_coeffs(:,i) = temp_coeffs(temp_coeffs>0);
    nu.grad_constraint(:,:,i) = temp_constr(any(temp_constr,2),:);  
end

%%

%Set up Hessian of gamma(x)

for i=1:n
   for j=1:n
        temp = gamma.grad_monoms(:,i) - ones(Coeff_total(n,gamma.deg-1),1);
        temp_coeffs = gamma.grad_monoms(:,i).*gamma.grad_coeffs(:,j);
        s = gamma.grad_monoms(:,i);
        for u=1:length(s)
            temp_constraint(u,:) = gamma.grad_constraint(u,:,j)*s(u);
        end
        t = temp_constraint;
        tt=t(any(t,2),:);
        gamma.hess_monoms(:,i) = temp(temp>=0);
        gamma.hess_coeffs(i,j,:) = temp_coeffs(temp_coeffs>0);
        gamma.hess_constraint(:,:,i,j) = tt;
        gamma.hess_monom_coeffs = Coeff_total(n,gamma.deg-2);
         
   end   
end


%Set up Hessian of nu(x)
for i=1:n
   for j=1:n
        temp = nu.grad_monoms(:,i) - ones(Coeff_total(n,nu.deg-1),1);
        temp_coeffs = nu.grad_monoms(:,i).*nu.grad_coeffs(:,j);
                s = nu.grad_monoms(:,i);
        for u=1:length(s)
            temp_constraint(u,:) = nu.grad_constraint(u,:,j)*s(u);
        end
        t = temp_constraint;
        tt=t(any(t,2),:);
        nu.hess_monoms(:,i) = temp(temp>=0);
        nu.hess_coeffs(i,j,:) = temp_coeffs(temp_coeffs>0);
        nu.hess_constraint(:,:,i,j) = tt;
        nu.hess_monom_coeffs = Coeff_total(n,gamma.deg-2);
         
   end   
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding Handelman representations for the quadratic forms of gamma(x) and
%nu(x).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Set up first polytope
% Gamma_1 = [0 0 1 ; 0 1 0 ; 1 0 -1 ; 1 -1 0];
% n_gamma_1 = 4;
% %Set up second polytope
% Gamma_2 = [0 0 1 ; 0 1 0 ; 1 0 -1 ; 1 -1 0];
% n_gamma_2 = 4;
% % %Set up third polytope (for LP solution only)
% % Gmmma_3 = [0 0 0 1 1; 0 0 1 1 0; 0 0 1 -1 0];

r = 1;
%Set up first polytope
Gamma_1 = [r 0 1 ; r 1 0 ; r 0 -1 ; r -1 0];
n_gamma_1 = 4;
%Set up second polytope
Gamma_2 = [r 0 1 ; r 1 0 ; r 0 -1 ; r -1 0];
n_gamma_2 = 4;
% %Set up third polytope (for LP solution only)
% Gmmma_3 = [0 0 0 1 1; 0 0 1 1 0; 0 0 1 -1 0];


%Set up Handelman form for Hessian of gamma(x)
disp('Set up Handelman form for Hessian of gamma(x)... ')
n_states = n;
V_gamma.monom_exps = lex_exps(n,gamma.deg-2);
V_gamma.poly_exps = lex_exps(n_gamma_1,gamma.deg-2);
V_gamma.poly_coeffs = Coeff_total(n_gamma_1,gamma.deg-2);
V_gamma.monom_coeffs = Coeff_total(n,gamma.deg-2);
V_gamma.temp_coeffs = zeros(V_gamma.monom_coeffs,V_gamma.poly_coeffs);

for i=1:V_gamma.poly_coeffs
    temp_cell = cell(n_gamma_1,1);
    tempexp = V_gamma.poly_exps(i,:);
    for j=1:n_gamma_1
        temp_cell{j} = multinomial(Gamma_1(j,:),tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    V_gamma.temp_coeffs(1:K,i) = temp_coeff;
end

%Set up Handelman form for Hessian of nu(x)
disp('Set up Handelman form for Hessian of nu(x)... ')
offset = 0;
V_nu.monom_exps = lex_exps(n,nu.deg-2+offset);
V_nu.poly_exps = lex_exps(n_gamma_2,nu.deg-2+offset);
V_nu.poly_coeffs = Coeff_total(n_gamma_2,nu.deg-2+offset);
V_nu.monom_coeffs = Coeff_total(n,nu.deg-2+offset);
V_nu.temp_coeffs = zeros(V_nu.monom_coeffs,V_nu.poly_coeffs);

for i=1:V_nu.poly_coeffs
    temp_cell = cell(n_gamma_2,1);
    tempexp = V_nu.poly_exps(i,:);
    for j=1:n_gamma_2
        temp_cell{j} = multinomial(Gamma_2(j,:),tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    V_nu.temp_coeffs(1:K,i) = temp_coeff;
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC Decomposition using Semidefinite Programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if((strcmp(relax_type,'SDP')==1))
    disp('Proceeding to set up and solve Handelman based Semidefinite Program...')
end


%
%Finding indices of nonzero coefficients of Hessian matrix.
disp('Finding indices of nonzero coefficients of Hessian matrix of gamma(x)...')


for k=1:gamma.hess_monom_coeffs 
    for l=1:gamma.coeff_total
        gamma.hess_matcoeffs(:,:,k,l) = zeros(n);  
        for i=1:n
            for j=1:n 
                gamma.hess_matcoeffs(i,j,k,l) = gamma.hess_constraint(k,l,i,j);
            end
        end
    end 
end





disp('Finding indices of nonzero coefficients of Hessian matrix of nu(x)...')

for k=1:nu.hess_monom_coeffs 
    for l=1:nu.coeff_total
        nu.hess_matcoeffs(:,:,k,l) = zeros(n);  
        for i=1:n
            for j=1:n 
                nu.hess_matcoeffs(i,j,k,l) = nu.hess_constraint(k,l,i,j);
            end
        end
    end 
end

%

%%
I = eye(P.coeff_total);
Z1 = zeros(P.coeff_total);
Z = zeros(6,15);
M = -V_nu.temp_coeffs;
Ag1m = gamma.hess_constraint(:,:,1,1)-gamma.hess_constraint(:,:,1,2); 
Ag1p = gamma.hess_constraint(:,:,1,1)+gamma.hess_constraint(:,:,1,2);
Ag2m = gamma.hess_constraint(:,:,2,2)-gamma.hess_constraint(:,:,2,1); 
Ag2p = gamma.hess_constraint(:,:,2,2)+gamma.hess_constraint(:,:,2,1);
Ah1m = nu.hess_constraint(:,:,1,1)-nu.hess_constraint(:,:,1,2); 
Ah1p = nu.hess_constraint(:,:,1,1)+nu.hess_constraint(:,:,1,2);
Ah2m = nu.hess_constraint(:,:,2,2)-nu.hess_constraint(:,:,2,1); 
Ah2p = nu.hess_constraint(:,:,2,2)+nu.hess_constraint(:,:,2,1);
nudiag=  nu.hess_constraint(:,:,1,1)+nu.hess_constraint(:,:,2,2);
rho = zeros(69);
rho(64) = 1;

 Aeq = [I,-I,Z1,Z1,Z1,Z1,Z1,Z1,Z1,Z1,Z1;Ag1m,Z,M,Z,Z,Z,Z,Z,Z,Z,Z;Ag1p,Z,Z,M,Z,Z,Z,Z,Z,Z,Z;Ag2m,Z,Z,Z,M,Z,Z,Z,Z,Z,Z;Ag2p,Z,Z,Z,Z,M,Z,Z,Z,Z,Z;Z,Ah1m,Z,Z,Z,Z,M,Z,Z,Z,Z;Z,Ah1p,Z,Z,Z,Z,Z,M,Z,Z,Z;Z,Ah2m,Z,Z,Z,Z,Z,Z,M,Z,Z;Z,Ah2p,Z,Z,Z,Z,Z,Z,Z,M,Z;Z,-nudiag,Z,Z,Z,Z,Z,Z,Z,Z,M];
[r,l] = size(Aeq)

rho = zeros(69,1);
rho(64) = 1;

 Aeq = [Aeq,rho];
 % Aeq = [I,-I,Z1,Z1,Z1,Z1,Z1,Z1,Z1,Z1;Ag1m,Z,M,Z,Z,Z,Z,Z,Z,Z;Ag1p,Z,Z,M,Z,Z,Z,Z,Z,Z;Ag2m,Z,Z,Z,M,Z,Z,Z,Z,Z;Ag2p,Z,Z,Z,Z,M,Z,Z,Z,Z;Z,Ah1m,Z,Z,Z,Z,M,Z,Z,Z;Z,Ah1p,Z,Z,Z,Z,Z,M,Z,Z;Z,Ah2m,Z,Z,Z,Z,Z,Z,M,Z;Z,Ah2p,Z,Z,Z,Z,Z,Z,Z,M];
[r,l] = size(Aeq);
beq = zeros(r,1);
beq(1:15) = P.coeffs;
f = zeros(l,1);
f(l) = 0;
Aineq = blkdiag(Z1,Z1,-I,-I,-I,-I,-I,-I,-I,-I,-I,-1);
bineq = zeros(166,1);
 bineq(31:150) = -0.000001*ones(120,1);
 [X,fval,exitflag] = linprog(f,Aineq,bineq,Aeq,beq);
% gcoeffs = XX(1:15);
% hcoeffs = XX(16:30);
% % XX = sedumi(Aeq,beq,c);
% 
% %  for i=1:gamma.coeff_total
% %      Gamma_Final_Coeffs(i) = value(gamma_sdpvars{i});
% %  end
% %  
% %   for i=1:nu.coeff_total
% %      Nu_Final_Coeffs(i) = value(nu_sdpvars{i});
% %   end
% %  
% %  
% %  test = Gamma_Final_Coeffs-Nu_Final_Coeffs;
% %  P.coeffs-test;
% %%
exitflag
fval
X1=X(1:15)
X2=X(16:30)
X1-X2

f(l) = 1;
Aineq = blkdiag(Z1,Z1,-I,-I,-I,-I,-I,-I,-I,-I,-I,-1);
bineq = zeros(166,1);
 bineq(31:150) = -0.000001*ones(120,1);
 [Y,fval,exitflag] = linprog(f,Aineq,bineq,Aeq,beq);
Y1=Y(1:15)
Y2=Y(16:30)
Y1-Y2
 
syms x y
for i=1:15
exp = P.monoms(i,:);
monoms(i) = (x^(exp(1)))*(y^exp(2));
end

f =  monoms*P.coeffs.';
g1 =  monoms*X1;
h1 =  monoms*X2;


g2 =  monoms*Y1;
h2 =  monoms*Y2;

figure

fsurf(f,[-1,1])

figure
hold on
fsurf(f,[-1,1])
fsurf(g1,[-1,1])
fsurf(-h1,[-1,1])

figure
hold on
fsurf(f,[-1,1])
fsurf(g2,[-1,1])
fsurf(-h2,[-1,1])
% ff = subs(f,y,0);
% gg = subs(g,y,0);
% hh = subs(h,y,0);
% % hold on
% figure 
% hold on
% fplot(ff,[-1,1]);
% fplot(gg,[-1,1]);
% fplot(hh,[-1,1]);
% 
% figure 
% hold on
% 
% fplot(gg,[-1,1]);
% fplot(hh,[-1,1]);
% fplot(gg-hh,[-1,1])

toc
