
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
n = 3;                                                  %Dimension of state space





P.deg = 6;                                              %Degree of polynomial to be decomposed
P.coeff_total = Coeff_total(n,P.deg);                       %Total number of coefficients of P(x)
offset_deg = 0;                                         %Degree offset

gamma.deg = P.deg+offset_deg;                           %Degree of gamma(x)
nu.deg = P.deg+offset_deg;                              %Degree of nu(x)


P.monoms = lex_exps(n,P.deg);                           %Monomials of P(x)
P.coeffs = ones(1,Coeff_total(n,P.deg));               %Set up a vector of coefficients for P(x)

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

gammma.grad_monoms = zeros(n,Coeff_total(n,gamma.deg-1));       %Set up monomials for gradient of gamma(x).
nu.grad_monoms = zeros(n,Coeff_total(n,nu.deg-1))';             %Set up  monomials for gradient of nu(x).

gammma.hess_monoms = zeros(n,Coeff_total(n,gamma.deg-2));       %Set up monomials for gradient of gamma(x).
nu.hess_monoms = zeros(n,Coeff_total(n,nu.deg-2))';             %Set up monomials for gradient of nu(x).

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

%Set up first polytope
Gamma_1 = [0 0 1 0 ; 0 1 0 0; 1 0 -1 0 ; 1 -1 0 0; 0 0 0 1; 1 0 0 -1];
n_gamma_1 = 6;
%Set up second polytope
Gamma_2 = [0 0 1 0 ; 0 1 0 0; 1 0 -1 0 ; 1 -1 0 0; 0 0 0 1; 1 0 0 -1];
n_gamma_2 = 6;
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

V_nu.monom_exps = lex_exps(n,nu.deg-2);
V_nu.poly_exps = lex_exps(n_gamma_2,nu.deg-2);
V_nu.poly_coeffs = Coeff_total(n_gamma_2,nu.deg-2);
V_nu.monom_coeffs = Coeff_total(n,nu.deg-2);
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


%%
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



% 

disp('Finding indices of nonzero coefficients of Hessian matrix of gamma(x)...')

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

%%
% 
%Solving the SDP with YALMIP

I=eye(n);
Epsilon= 0.000001;
F = [];
for i=1:gamma.coeff_total
    gamma_sdpvars{i} = sdpvar(1,1);
end

for i=1:nu.coeff_total
    nu_sdpvars{i} = sdpvar(1,1);
end

for i=1:V_gamma.poly_coeffs
    gamma_hess_sdpvars{i} = sdpvar(n,n);
    F = [F,gamma_hess_sdpvars{i} - Epsilon*I>=0];
end

for i=1:V_nu.poly_coeffs
    nu_hess_sdpvars{i} = sdpvar(n,n);
    F = [F,nu_hess_sdpvars{i} - Epsilon*I>=0];
end

for i=1:P.coeff_total
    F = [F,gamma_sdpvars{i}-nu_sdpvars{i} == P.coeffs(i)];
end

for j=1:gamma.hess_monom_coeffs
    Const_gamma{j} = 0;
    for i=1:gamma.coeff_total
        Const_gamma{j} = Const_gamma{j}+gamma.hess_matcoeffs(:,:,j,i)*gamma_sdpvars{i};
    end
end

for j=1:nu.hess_monom_coeffs
    Const_nu{j} = 0;
    for i=1:nu.coeff_total
        Const_nu{j} = Const_nu{j}+nu.hess_matcoeffs(:,:,j,i)*nu_sdpvars{i};
    end
end

for j=1:V_gamma.monom_coeffs
    Const_handel_gamma{j} = 0;
    for i=1:V_gamma.poly_coeffs
        Const_handel_gamma{j} = Const_handel_gamma{j}+V_gamma.temp_coeffs(j,i)*gamma_hess_sdpvars{i};
    end
end

for j=1:V_nu.monom_coeffs
    Const_handel_nu{j} = 0;
    for i=1:V_nu.poly_coeffs
        Const_handel_nu{j} = Const_handel_nu{j}+V_nu.temp_coeffs(j,i)*nu_hess_sdpvars{i};
    end
end

for j=1:V_gamma.monom_coeffs
    F = [F, Const_gamma{j} == Const_handel_gamma{j}];
end

for j=1:V_nu.monom_coeffs
    F = [F, Const_nu{j} == Const_handel_nu{j}];
end


% F = [F, Const_gamma{1}== Const_handel_gamma{1}];
% F = [F, Const_gamma{2}== Const_handel_gamma{2}];
% F = [F, Const_gamma{3}==Const_handel_gamma{3}];
% 
% F = [F, Const_nu{1}== Const_handel_nu{1}];
% F = [F, Const_nu{2}== Const_handel_nu{2}];
% F = [F, Const_nu{3}==Const_handel_nu{3}];


% F = [F, gamma.hess_matcoeffs(:,:,1,4)*gamma_sdpvars{4} + gamma.hess_matcoeffs(:,:,1,5)*gamma_sdpvars{5} + gamma.hess_matcoeffs(:,:,1,6)*gamma_sdpvars{6}== V_gamma.temp_coeffs(1,1)*gamma_hess_sdpvars{1} + V_gamma.temp_coeffs(1,2)*gamma_hess_sdpvars{2} + V_gamma.temp_coeffs(1,3)*gamma_hess_sdpvars{3}];
% F = [F, gamma.hess_matcoeffs(:,:,2,7)*gamma_sdpvars{7} + gamma.hess_matcoeffs(:,:,2,8)*gamma_sdpvars{8}+gamma.hess_matcoeffs(:,:,2,9)*gamma_sdpvars{9}== V_gamma.temp_coeffs(2,3)*gamma_hess_sdpvars{3} + V_gamma.temp_coeffs(2,5)*gamma_hess_sdpvars{5}];
% F = [F, gamma.hess_matcoeffs(:,:,3,8)*gamma_sdpvars{8} + gamma.hess_matcoeffs(:,:,3,9)*gamma_sdpvars{9} + gamma.hess_matcoeffs(:,:,3,10)*gamma_sdpvars{10}==V_gamma.temp_coeffs(3,2)*gamma_hess_sdpvars{2} + V_gamma.temp_coeffs(3,4)*gamma_hess_sdpvars{4}];

% F = [F, nu.hess_matcoeffs(:,:,1,4)*nu_sdpvars{4} + nu.hess_matcoeffs(:,:,1,5)*nu_sdpvars{5} + nu.hess_matcoeffs(:,:,1,6)*nu_sdpvars{6}== V_nu.temp_coeffs(1,1)*nu_hess_sdpvars{1} + V_nu.temp_coeffs(1,2)*nu_hess_sdpvars{2} + V_nu.temp_coeffs(1,3)*nu_hess_sdpvars{3}];
% F = [F, nu.hess_matcoeffs(:,:,2,7)*nu_sdpvars{7} + nu.hess_matcoeffs(:,:,2,8)*nu_sdpvars{8}+nu.hess_matcoeffs(:,:,2,9)*nu_sdpvars{9}== V_nu.temp_coeffs(2,3)*nu_hess_sdpvars{3} + V_nu.temp_coeffs(2,5)*nu_hess_sdpvars{5}];
% F = [F, nu.hess_matcoeffs(:,:,3,8)*nu_sdpvars{8} + nu.hess_matcoeffs(:,:,3,9)*nu_sdpvars{9} + nu.hess_matcoeffs(:,:,3,10)*nu_sdpvars{10}==V_nu.temp_coeffs(3,2)*nu_hess_sdpvars{2} + V_nu.temp_coeffs(3,4)*nu_hess_sdpvars{4}];


 optimize(F);
 
 for i=1:gamma.coeff_total
     Gamma_Final_Coeffs(i) = value(gamma_sdpvars{i});
 end
 
  for i=1:nu.coeff_total
     Nu_Final_Coeffs(i) = value(nu_sdpvars{i});
 end
 
%%
 for i=1:V_gamma.poly_coeffs
     Gamma_Handel_Final_Coeffs{i} = value(gamma_hess_sdpvars{i});
     A{i} = eigs(Gamma_Handel_Final_Coeffs{i});
 end
 for i=1:V_nu.poly_coeffs
     Nu_Handel_Final_Coeffs{i} = value(nu_hess_sdpvars{i});
     B{i} = eigs(Gamma_Handel_Final_Coeffs{i});
 end
toc


