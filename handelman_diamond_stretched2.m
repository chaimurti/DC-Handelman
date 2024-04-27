clc 
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stability Analysis on polytopes using piecewise polynomial Lyapunov
%functions synthesized using Handelman forms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
disp('Setting up the polytope and Lyapunov functions');

disp('Initial setup')

V.deg = 8;                                          %Degree of Lyapunov function
V1.deg = V.deg;                                          %Degree of Lyapunov function
V2.deg = V.deg;                                          %Degree of Lyapunov function
V3.deg = V.deg;                                          %Degree of Lyapunov function
V4.deg = V.deg;                                          %Degree of Lyapunov function
A.deg = 3;                                          %Degree of Vec field
Vdot.deg = V.deg+A.deg-1;                            %Degree of Vdot
V1dot.deg = V.deg+A.deg-1;                            %Degree of Vdot
V2dot.deg = V.deg+A.deg-1;                            %Degree of Vdot
V3dot.deg = V.deg+A.deg-1;                            %Degree of Vdot
V4dot.deg = V.deg+A.deg-1;                            %Degree of Vdot
n_polytope = 3;                                     %Number of sides in the polytope
n_states = 2;                                       %No. of state variables
V.monom_coeffs = Coeff_total(V.deg,n_states);       %Coeffs in monom. basis
V.poly_coeffs = Coeff_total(V.deg,n_polytope);      %Coeffs in poly. basis
Vdot.monom_coeffs = Coeff_total(Vdot.deg,n_states);
Vdot.poly_coeffs = Coeff_total(Vdot.deg,n_polytope);
V1.monom_coeffs = V.monom_coeffs;
V2.monom_coeffs = V.monom_coeffs;
V3.monom_coeffs = V.monom_coeffs;
V4.monom_coeffs = V.monom_coeffs;
V1.poly_coeffs = V.poly_coeffs;
V2.poly_coeffs = V.poly_coeffs;
V3.poly_coeffs = V.poly_coeffs;
V4.poly_coeffs = V.poly_coeffs;
V1dot.monom_coeffs = Vdot.monom_coeffs ;
V2dot.monom_coeffs = Vdot.monom_coeffs ;
V3dot.monom_coeffs = Vdot.monom_coeffs ;
V4dot.monom_coeffs = Vdot.monom_coeffs ;
V1dot.poly_coeffs = Vdot.poly_coeffs ;
V2dot.poly_coeffs = Vdot.poly_coeffs ;
V3dot.poly_coeffs = Vdot.poly_coeffs ;
V4dot.poly_coeffs = Vdot.poly_coeffs ;

epsilon = 0.01;

% %Set up first polytope
% l1.v{1} = [1 -1 0];
% l1.v{2} = [0 1 1];
% l1.v{3} = [0 1 -1];
% 
% %Set up second polytope
% l2.v{1} = [1 0 -1];
% l2.v{2} = [0 -1 1];
% l2.v{3} = [0 1 1];
% 
% %Set up third polytope
% l3.v{1} = [1 1 0];
% l3.v{2} = [0 -1 1];
% l3.v{3} = [0 -1 -1];
% 
% %Set up fourth polytope
% l4.v{2} = [0 1 -1];
% l4.v{1} = [1 0 1];
% l4.v{3} = [0 -1 -1];


% r1 = 2.25;   % worked only with linprog, deg=6
% r2 = 2.25;

% r1 = 2.4;   % worked only with linprog, deg=8
% r2 = 2.4;

% r1 = 2.3;   % worked only with linprog & SeDuMi, deg=8
% r2 = 2.5;

r1 = 2.5;   
r2 = 2.9;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set up first polytope
l1.v{1} = [r1*r2 -r2 -r1];
l1.v{2} = [0 1 0];
l1.v{3} = [0 0 1];

%Set up second polytope
l2.v{1} = [r1*r2 -r2 r1];
l2.v{2} = [0 1 0];
l2.v{3} = [0 0 -1];

%Set up third polytope
l3.v{1} = [r1*r2 r2 r1];
l3.v{2} = [0 -1 0];
l3.v{3} = [0 0 -1];

%Set up fourth polytope
l4.v{1} = [r1*r2 r2 -r1];
l4.v{2} = [0 -1 0];
l4.v{3} = [0 0 1];

%%
%Set up the vector field to be analyzed:
disp('Setting up f(x)')


A.exp = lex_exps(n_states,A.deg);
A.coeff_total = Coeff_total(n_states,A.deg);
A.coeffs = zeros(n_states,Coeff_total(n_states,A.deg));

k = lex_index_nh([0 1]);
A.coeffs(1,k) = -1;
A.coeffs(2,k) = -1;
j = lex_index_nh([0 1]);
k = lex_index_nh([1 0]);
A.coeffs(j,k) = 1;
k = lex_index_nh([2 1]);
A.coeffs(j,k) = 1;
% 

%%
%Testing which facets have nonzero coeffs
facet_test_vec = [1 0 0].';
facet_rule_1 = zeros(size(l1.v));
for i=1:length(l1.v)
    I = l1.v{i}*facet_test_vec;
    if I==0
        facet_rule_1(i) = 1;
    end
end

facet_rule_2 = zeros(size(l2.v));
for i=1:length(l2.v)
    I = l2.v{i}*facet_test_vec;
    if I==0
        facet_rule_2(i) = 1;
    end
end

facet_rule_3 = zeros(size(l3.v));
for i=1:length(l3.v)
    I = l3.v{i}*facet_test_vec;
    if I==0
        facet_rule_3(i) = 1;
    end
end

facet_rule_4 = zeros(size(l4.v));
for i=1:length(l4.v)
    I = l4.v{i}*facet_test_vec;
    if I==0
        facet_rule_4(i) = 1;
    end
end



%%
%Setting up V1 in terms of g(b_alpha)

disp('Setting up V1 in terms of g(b_alpha)')
 
V1.monom_exps = lex_exps(n_states,V.deg);
V1.poly_exps = lex_exps(n_polytope, V.deg);
G1.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V1.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l1.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G1.temp_coeffs(1:K,i) = temp_coeff;
end

[G1.coeffs,G1.exps,G1.coeff_total] = coeff_elim(G1.temp_coeffs,V1.poly_exps,facet_rule_1);

%%

%Setting up V2 in terms of g(b_alpha)

disp('Setting up V2 in terms of g(b_alpha)')
 
V2.monom_exps = lex_exps(n_states,V.deg);
V2.poly_exps = lex_exps(n_polytope, V.deg);
G2.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V2.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V2.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l2.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G2.temp_coeffs(1:K,i) = temp_coeff;
end

[G2.coeffs,G2.exps,G2.coeff_total] = coeff_elim(G2.temp_coeffs,V2.poly_exps,facet_rule_2);

%%
%Setting up V3 in terms of g(b_alpha)

disp('Setting up V3 in terms of g(b_alpha)')
 
V3.monom_exps = lex_exps(n_states,V.deg);
V3.poly_exps = lex_exps(n_polytope, V.deg);
G3.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V3.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V3.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l3.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G3.temp_coeffs(1:K,i) = temp_coeff;
end

[G3.coeffs,G3.exps,G3.coeff_total] = coeff_elim(G3.temp_coeffs,V3.poly_exps,facet_rule_3);

%%

%Setting up V4 in terms of g(b_alpha)

disp('Setting up V4 in terms of g(b_alpha)')
 
V4.monom_exps = lex_exps(n_states,V.deg);
V4.poly_exps = lex_exps(n_polytope, V.deg);
G4.temp_coeffs = zeros(V.monom_coeffs,V.poly_coeffs);
for i=1:V4.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V4.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l4.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    G4.temp_coeffs(1:K,i) = temp_coeff;
end

[G4.coeffs,G4.exps,G4.coeff_total] = coeff_elim(G4.temp_coeffs,V4.poly_exps,facet_rule_4);


%%
%Setting up V1dot in terms of l(c_alpha)

disp('Setting up V1dot in terms of l(c_alpha)')

V1dot.monom_exps = lex_exps(n_states,Vdot.deg);
V1dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L1.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V1dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V1dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l1.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L1.coeffs(1:K,i) = -temp_coeff;
end

[L1.coeffs,L1.exps,L1.coeff_total] = coeff_elim(L1.coeffs,V1dot.poly_exps,facet_rule_1);


%%
%Setting up V2dot in terms of l(c_alpha)

disp('Setting up V2dot in terms of l(c_alpha)')

V2dot.monom_exps = lex_exps(n_states,Vdot.deg);
V2dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L2.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V2dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V2dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l2.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L2.coeffs(1:K,i) = -temp_coeff;
end

[L2.coeffs,L2.exps,L2.coeff_total] = coeff_elim(L2.coeffs,V2dot.poly_exps,facet_rule_2);


%%
%Setting up V3dot in terms of l(c_alpha)

disp('Setting up V3dot in terms of l(c_alpha)')

V3dot.monom_exps = lex_exps(n_states,Vdot.deg);
V3dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L3.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V3dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V3dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l3.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L3.coeffs(1:K,i) = -temp_coeff;
end

[L3.coeffs,L3.exps,L3.coeff_total] = coeff_elim(L3.coeffs,V3dot.poly_exps,facet_rule_3);


%%
%Setting up V4dot in terms of l(c_alpha)

disp('Setting up V4dot in terms of l(c_alpha)')

V4dot.monom_exps = lex_exps(n_states,Vdot.deg);
V4dot.poly_exps = lex_exps(n_polytope, Vdot.deg);
L4.coeffs = zeros(Vdot.monom_coeffs,Vdot.poly_coeffs);
for i=1:V4dot.poly_coeffs
    temp_cell = cell(n_polytope,1);
    tempexp = V4dot.poly_exps(i,:);
    for j=1:n_polytope
        temp_cell{j} = multinomial(l4.v{j},tempexp(j));
    end
    temp_coeff = polyprod2(temp_cell,n_states,tempexp);
    K = length(temp_coeff);
    L4.coeffs(1:K,i) = -temp_coeff;
end

[L4.coeffs,L4.exps,L4.coeff_total] = coeff_elim(L4.coeffs,V4dot.poly_exps,facet_rule_4);

Vdot.poly_coeffs = L4.coeff_total;

%%
%                                                                                                                                                                                                                                                                                                                               
H1.coeffs = zeros(V1dot.monom_coeffs,G1.coeff_total);
GradV1.exps = lex_exps(n_states,V.deg-1);
GradV1.coeff_total = Coeff_total(n_states,V.deg-1);

exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV1.coeff_total,G1.coeff_total);
    for j=1:V1.monom_coeffs
        temp_exp_1 = V1.monom_exps(j,:);
        temp_coeff_vec = G1.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k=lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    temp_coeffs;
    I = 4-lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
%       pause
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        temp_A = A_temp(j);
        for l=1:GradV1.coeff_total
            temp_exp_4 = temp_exp_3+GradV1.exps(l,:);
            temp_coeff_h = temp_A*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H1.coeffs(k,:) = H1.coeffs(k,:)+temp_coeff_h;
        end
    end
end

%%
%                                                                                                                                                                                                                                                                                                                               
H2.coeffs = zeros(Vdot.monom_coeffs,G2.coeff_total);
GradV2.exps = lex_exps(n_states,V.deg-1);
GradV2.coeff_total = Coeff_total(n_states,V.deg-1);
GradV2.prodcoeffs = ones(2,GradV2.coeff_total);
exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV2.coeff_total,G2.coeff_total);
    for j=1:V2.monom_coeffs
        temp_exp_1 = V2.monom_exps(j,:);
        temp_coeff_vec = G2.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k=lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    
    I = 4 -lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
    
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        for l=1:GradV2.coeff_total
            temp_exp_4 = temp_exp_3+GradV2.exps(l,:);
            temp_coeff_h = A_temp(j)*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H2.coeffs(k,:) = H2.coeffs(k,:)+temp_coeff_h;
        end
    end
end


%%
%                                                                                                                                                                                                                                                                                                                               
H3.coeffs = zeros(V3dot.monom_coeffs,G3.coeff_total);
GradV3.exps = lex_exps(n_states,V.deg-1);
GradV3.coeff_total = Coeff_total(n_states,V.deg-1);
GradV3.prodcoeffs = ones(2,GradV3.coeff_total);
exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV3.coeff_total,G3.coeff_total);
    for j=1:V3.monom_coeffs
        temp_exp_1 = V3.monom_exps(j,:);
        temp_coeff_vec = G3.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k=lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    
    I = 4-lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
    
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        for l=1:GradV3.coeff_total
            temp_exp_4 = temp_exp_3+GradV3.exps(l,:);
            temp_coeff_h = A_temp(j)*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H3.coeffs(k,:) = H3.coeffs(k,:)+temp_coeff_h;
        end
    end
end

%%
%                                                                                                                                                                                                                                                                                                                               
H4.coeffs = zeros(Vdot.monom_coeffs,G4.coeff_total);
GradV4.exps = lex_exps(n_states,V.deg-1);
GradV4.coeff_total = Coeff_total(n_states,V.deg-1);
GradV4.prodcoeffs = ones(2,GradV4.coeff_total);
exp_slider_mat = eye(n_states);


for i=1:n_states
    exp_slider = exp_slider_mat(i,:);
    temp_coeffs = zeros(GradV4.coeff_total,G4.coeff_total);
    for j=1:V4.monom_coeffs
        temp_exp_1 = V4.monom_exps(j,:);
        temp_coeff_vec = G4.coeffs(j,:);
        if temp_exp_1(i)~=0
            temp_exp_2 = temp_exp_1-exp_slider;
            temp_coeff_vec_2 = temp_coeff_vec*temp_exp_1(i);
            k=lex_index_nh(temp_exp_2);
            temp_coeffs(k,:) = temp_coeff_vec_2;
        end
    end
    
    I = 4-lex_index_nh(exp_slider);
    A_temp = A.coeffs(I,:);
    
    for j=1:A.coeff_total
        temp_exp_3 = A.exp(j,:);
        for l=1:GradV4.coeff_total
            temp_exp_4 = temp_exp_3+GradV4.exps(l,:);
            temp_coeff_h = A_temp(j)*temp_coeffs(l,:);
            k= lex_index_nh(temp_exp_4);
            H4.coeffs(k,:) = H4.coeffs(k,:)+temp_coeff_h;
        end
    end
end


%%
%Setting up blocks for interface between P1 and P2

G12.coeffs = coeff_elim_zeros(G1.coeffs,G1.exps,[0 0 1]);
G21.coeffs = -coeff_elim_zeros(G2.coeffs,G2.exps,[0 0 1]);

%%
%Setting up blocks for interface between P2 and P3

G23.coeffs = coeff_elim_zeros(G2.coeffs,G2.exps,[0 1 0]);
G32.coeffs = -coeff_elim_zeros(G3.coeffs,G3.exps,[0 1 0]);

%%
%Setting up blocks for interface between P3 and P4

G34.coeffs = coeff_elim_zeros(G3.coeffs,G3.exps,[0 0 1]);
G43.coeffs = -coeff_elim_zeros(G4.coeffs,G4.exps,[0 0 1]);

%%
%Setting up blocks for interface between P4 and P1

G41.coeffs = coeff_elim_zeros(G4.coeffs,G4.exps,[0 1 0]);
G14.coeffs = -coeff_elim_zeros(G1.coeffs,G1.exps,[0 1 0]);

%%
%  %Set up Matrices for LP
% % 
% % 
% 

 disp('Setting up and solving the LP')
% 
% Vd1 = [H1.coeffs,-L1.coeffs];
% Vd2 = [H2.coeffs,-L2.coeffs];
% Vd3 = [H3.coeffs,-L3.coeffs];
% Vd4 = [H4.coeffs,-L4.coeffs];
% Amat_temp = zeros(0,0);
% Amat_temp = blkdiag(Vd1,Vd2,Vd3,Vd4);
% Ll = size(Amat_temp);
%  block1 = zeros(V.monom_coeffs,Vdot.poly_coeffs);
%  block2 = zeros(V.monom_coeffs,Vdot.poly_coeffs+G1.coeff_total);
%  Cont12 = [G12.coeffs,block1,G21.coeffs,block1,block2,block2];
%  Cont23 = [block2,G23.coeffs,block1,G32.coeffs,block1,block2];
%  Cont34 = [block2,block2,G34.coeffs,block1,G43.coeffs,block1];
%  Cont41 = [G14.coeffs,block1,block2,block2,G41.coeffs,block1];
%  Amat_temp = [Amat_temp;Cont12;Cont23;Cont34;Cont41];
%  L2 = size(Amat_temp);
%  Block1 = [eye(G1.coeff_total);zeros(3*G1.coeff_total,G1.coeff_total)];
%  Block2 = [zeros(G1.coeff_total);eye(G1.coeff_total);zeros(2*G1.coeff_total,G1.coeff_total)];
%  Block3 = [zeros(2*G1.coeff_total,G1.coeff_total);eye(G1.coeff_total);zeros(G1.coeff_total)];
%  Block4 = [zeros(3*G1.coeff_total,G1.coeff_total);eye(G1.coeff_total)];
%  BlockC = zeros(4*G1.coeff_total,Vdot.poly_coeffs);
%  BlockA = zeros(L2(1),4*G1.coeff_total);
%  BlockB = [Block1,BlockC,Block2,BlockC,Block3,BlockC,Block4,BlockC];
% BlockD= [zeros(L2(1),4*G1.coeff_total);-eye(4*G1.coeff_total)];
% Amat_temp = [Amat_temp;BlockB];
% Amat = [Amat_temp,BlockD];
% L3 = size(Amat);
% cmat = zeros(L3(2),1).';
% bmat = zeros(L3(1),1).';
% 
% for j=(L3(1)-4*G1.coeff_total+1):L3(1)
%     bmat(j) = epsilon;
% end

% 
% B1 = blkdiag(H1.coeffs,H2.coeffs,H3.coeffs,H4.coeffs);
% B2 = blkdiag(-L1.coeffs,-L2.coeffs,-L3.coeffs,-L4.coeffs);
% B3 = eye(4*G1.coeff_total);
% zeroblock = zeros(V.monom_coeffs,G1.coeff_total);
% B4 = [G12.coeffs,G21.coeffs,zeroblock,zeroblock;zeroblock,G23.coeffs,G32.coeffs,zeroblock;zeroblock,zeroblock,G34.coeffs,G43.coeffs;G14.coeffs,zeroblock,zeroblock,G41.coeffs];
% 
% A1 = [B1;B4;B3];
% A2 = [B2;zeros(4*V.monom_coeffs,4*Vdot.poly_coeffs);zeros(4*G1.coeff_total,4*Vdot.poly_coeffs)];
% A3 = [zeros(4*Vdot.monom_coeffs+4*V.monom_coeffs,4*G1.coeff_total);-B3];
% Amat = [A1,A2,A3];
% A1 = [B1,B2];
% A2 = [B4,zeros(4*V.monom_coeffs,Vdot.poly_coeffs)];
% A3 = [B3, zeros(4*G1.coeff_total,Vdot.poly_coeffs)];
% A = [A1;A2;A3]; 4*vdot.monom_coeffs+4*G1.coeff_total,G1.coeff_total);-B3]];


B1 = blkdiag(H1.coeffs,H2.coeffs,H3.coeffs,H4.coeffs);

B2 = blkdiag(-L1.coeffs,-L2.coeffs,-L3.coeffs,-L4.coeffs);

zeroblock = zeros(V.monom_coeffs,G1.coeff_total);
B41 = [G12.coeffs,G21.coeffs];
B42 = [G23.coeffs,G32.coeffs];
B43 = [G34.coeffs,G43.coeffs];
B44 = [G14.coeffs,zeros(V.monom_coeffs,2*G1.coeff_total),G41.coeffs];
B4 = [B41,zeros(V.monom_coeffs,2*G1.coeff_total);zeros(V.monom_coeffs,G1.coeff_total),B42,zeros(V.monom_coeffs,G1.coeff_total);zeros(V.monom_coeffs,2*G1.coeff_total),B43;B44];

A1 = [B1;B4];

A2 = [B2;zeros(4*V.monom_coeffs,4*Vdot.poly_coeffs)];

% A2 = [B2;zeros(V1.monom_coeffs,4*Vdot.poly_coeffs)];
% A2 = B2;
Aeq = [A1,A2];
Aeq = constraint_elim(Aeq);
A_ineq = eye(4*G1.coeff_total+4*Vdot.poly_coeffs); 

Aeq_size = size(Aeq);
cmat = zeros(Aeq_size(2),1);
bmat = zeros(Aeq_size(1),1);
b_ineq = zeros(4*G1.coeff_total+4*Vdot.poly_coeffs,1);

%Epsilon that worked for r<1 := (3,4,5)
%Epsilon that worked for r=1 := (3,5)

for j=1:4
    b_ineq(3+G1.coeff_total*(j-1)) = epsilon;
    b_ineq(4+G1.coeff_total*(j-1)) = 0;
    b_ineq(5+G1.coeff_total*(j-1)) = epsilon;
end

% options.Display='iter';
% [XX,fval,exitflag] = linprog(cmat,-A_ineq,-b_ineq,Aeq,bmat,[],[],[],options);

%options = optimset('LargeScale','off','Simplex','off');
% options.MaxIter=1000;
% options.Display='iter';
% [XX,fval,exitflag,output] = linprog(cmat,-A_ineq,-b_ineq,Aeq,bmat,[],[],[],options);

%[XX,y,info] = sedumi(Aeq,bmat,cmat);

solver_g
XX=result.x;

%   temp = 0;
%   X1 = X(temp+1:temp+G1.coeff_total);
%   temp = temp+ G1.coeff_total;
%   X2 = X(temp+1:G2.coeff_total+temp);
%   temp = temp+G2.coeff_total;
%   X3 = X(temp+1:G3.coeff_total+temp);
%   temp = temp+G3.coeff_total;
%   X4 = X(temp+1:G4.coeff_total+temp);
%   temp = temp+G4.coeff_total;
%   
%   C1=X(temp+1:Vdot.poly_coeffs+temp);
%   temp = temp+Vdot.poly_coeffs;
%     C2=X(temp+1:Vdot.poly_coeffs+temp);
%   temp = temp+Vdot.poly_coeffs;
%     C3=X(temp+1:Vdot.poly_coeffs+temp);
%   temp = temp+Vdot.poly_coeffs;
%     C4=X(temp+1:Vdot.poly_coeffs+temp);
%   temp = temp+Vdot.poly_coeffs;
%   

  disp ('Verification');
  temp = 0;
  XX1 = G1.coeffs*XX(temp+1:temp+G1.coeff_total);
  temp = temp+ G1.coeff_total;
  XX2 = G2.coeffs*XX(temp+1:G2.coeff_total+temp);
  temp = temp+G2.coeff_total;
  XX3 = G3.coeffs*XX(temp+1:G3.coeff_total+temp);
  temp = temp+G3.coeff_total;
  XX4 = G4.coeffs*XX(temp+1:G4.coeff_total+temp);
  temp = temp+G4.coeff_total;
  exps = lex_exps(2,V.deg);
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(XX1,2,exps,[i,j]);
          if test<0
              disp('Negative value');
              disp(test);
          end
      end
  end
  
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(XX2,2,exps,[i,-j]);
          if test<0
              disp('Negative value');
              disp(test);
          end
      end
  end
  
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(XX3,2,exps,[-i,-j]);
          if test<0
              disp('Negative value');
              disp(test);
          end
      end
  end
  
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(XX4,2,exps,[-i,j]);
          if test<0
              disp('Negative value');
              disp(test);
          end
      end
  end
  
  temp = 0;
  DX1 = H1.coeffs*XX(temp+1:temp+G1.coeff_total);
  temp = temp+ G1.coeff_total;
  DX2 = H2.coeffs*XX(temp+1:G2.coeff_total+temp);
  temp = temp+G2.coeff_total;
  DX3 = H3.coeffs*XX(temp+1:G3.coeff_total+temp);
  temp = temp+G3.coeff_total;
  DX4 = H4.coeffs*XX(temp+1:G4.coeff_total+temp);
  temp = temp+G4.coeff_total;
  exps = lex_exps(2,Vdot.deg);
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(DX1,2,exps,[i,j]);
          if test>0
              disp('Pos value');
              disp(test);
          end
      end
  end
  
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(DX2,2,exps,[i,-j]);
          if test>0
              disp('Pos value');
              disp(test);
          end
      end
  end
  
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(DX3,2,exps,[-i,-j]);
          if test>0
              disp('Pos value');
              disp(test);
          end
      end
  end
  
  for i=0:.01:r1
      for j=0:.01:(r2-(r2/r1)*i)
          test = lyap_verify2(DX4,2,exps,[-i,j]);
          if test>0
              disp('Pos value');
              disp(test);
          end
      end
  end
  
  temp = 0;
  GX12 = G12.coeffs*XX(temp+1:temp+G1.coeff_total);
  GX14 = G14.coeffs*XX(temp+1:temp+G1.coeff_total);
  temp = temp+ G1.coeff_total;
  GX21 = G21.coeffs*XX(temp+1:G2.coeff_total+temp);
  GX23 = G23.coeffs*XX(temp+1:G2.coeff_total+temp);
  temp = temp+G2.coeff_total;
  GX32 = G32.coeffs*XX(temp+1:G3.coeff_total+temp);
  GX34 = G34.coeffs*XX(temp+1:G3.coeff_total+temp);
  temp = temp+G3.coeff_total;
  GX43 = G43.coeffs*XX(temp+1:G4.coeff_total+temp);
  GX41 = G41.coeffs*XX(temp+1:G4.coeff_total+temp);
  temp = temp+G4.coeff_total;

  
  G12test = GX12+GX21;
  G23test = GX23+GX32;
  G34test = GX34+GX43;
  G41test = GX14+GX41;