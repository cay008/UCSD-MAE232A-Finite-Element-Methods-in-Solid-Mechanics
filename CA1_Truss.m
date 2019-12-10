close all; clear;clc

% Geometry Parameters
L = 10; H = 15;

A = [1 1 2 1.5 1.5 2 2]; % Cross Section of Areas
E = repmat(30e6,1,7); % Young's Modulus
Truss_L = [2*L 2*L norm([L,H]) norm([L,H]) norm([L,H]) norm([L,H]) 2*L]; % Truss Length
k_coeff = A.*E./Truss_L; % Local stiffness matrix

theta = [atan2(0,2*L) atan2(0,2*L) atan2(H,L) atan2(H,-L) atan2(H,L) atan2(H,-L) atan2(0,2*L)];

elements = [1 2; 2 3 ; 1 4; 2 4 ; 2 5 ; 3 5 ; 4 5];
nof_elements = size(elements,1);
nof_nodes = 5;

% Prelocations
R = zeros(2,2);
K = zeros(2*nof_nodes, 2*nof_nodes);
K_element = zeros(2*nof_nodes, 2*nof_nodes);
u = zeros(2*nof_nodes,1);
F = zeros(2*nof_nodes,1);

% Loading Condition
F_lc = [-20e3 ; -7.5e3]; % External Loading
F_lc_hz = [4]; F_lc_hz = 2*F_lc_hz -1 ;
F_lc_vt = [2]; F_lc_vt = 2*F_lc_vt;
F_lc_index = [F_lc_hz ; F_lc_vt];
F_lc_index = sort(F_lc_index,'ascend');
F(F_lc_index) = F_lc;

% Boundary Conditions
u_bc = [0 ; 0 ; 0];
u_bc_hz = [3]; u_bc_hz = 2*u_bc_hz -1 ;
u_bc_vt = [1 ; 3]; u_bc_vt = 2*u_bc_vt ;
u_bc_index = [u_bc_hz ; u_bc_vt];
u_bc_index = sort(u_bc_index,'ascend');
u(u_bc_index) = u_bc;

% Construciton of K (Global)
for i = 1 : nof_elements
    R = [cos(theta(i)).^2 cos(theta(i))*sin(theta(i)); cos(theta(i))*sin(theta(i)) sin(theta(i)).^2];
    K_element(2*elements(i,1) - 1 : 2*elements(i,1) , 2*elements(i,1) - 1 : 2*elements(i,1)) = k_coeff(i)*R;
    K_element(2*elements(i,1) - 1 : 2*elements(i,1) , 2*elements(i,2) - 1 : 2*elements(i,2)) = k_coeff(i)*(-R);
    K_element(2*elements(i,2) - 1 : 2*elements(i,2) , 2*elements(i,1) - 1 : 2*elements(i,1)) = k_coeff(i)*(-R);
    K_element(2*elements(i,2) - 1 : 2*elements(i,2) , 2*elements(i,2) - 1 : 2*elements(i,2)) = k_coeff(i)*R;
    K = K + K_element;
    K_element(:,:) = 0;
end
%% Q1. Nodal Displacements
disp('1. Displacements at All Nodes')
node_index_all = transpose(1 : 2*nof_nodes);
F_known_index = setdiff(node_index_all,u_bc_index,'stable');
u(F_known_index) = K(F_known_index,F_known_index) \ F(F_known_index)

% Reactions Forces
F = K*u;

%% Q2 Reactions Forces at nodes 1 and 3
nodes_index = [1,3];
nodes_index = [2*nodes_index(1) - 1 : 2*nodes_index(1), 2*nodes_index(2) - 1 : 2*nodes_index(2)];
F_13 = F(nodes_index);
disp('2. Reaction Forces at Node 1 and Node 2')
sprintf('Horizontal Reaction Force at Node 1 : %f lbf \n Vertical Reaction Force at Node 1 : %f lbf \n', F_13(1),F_13(2))
sprintf('Horizontal Reaction Force at Node 3 : %f lbf \n Vertical Reaction Force at Node 3 : %f lbf \n', F_13(3),F_13(4))

%% Q3. Stressess
axial_disp_2 = u(2*elements(:,2)-1).*cos(theta).' + u(2*elements(:,2)).*sin(theta).';
axial_disp_1 = u(2*elements(:,1)-1).*cos(theta).' + u(2*elements(:,1)).*sin(theta).';
axial_disp = axial_disp_2 - axial_disp_1;
strain = axial_disp./Truss_L.';
disp('3. Stress in elements (psi)')
stress = E.'.*strain


%% Q4. Optimizations
A(3:6) = 0;
elements_switch = perms([2 1.5 1.5 2]);
node2disp = zeros(size(elements_switch,1),1);
for j = 1 : size(elements_switch,1)
    A(3:6) = elements_switch(j,:);
    F_lc_hz = [4]; F_lc_vt = [2]; u_bc_hz = [3]; u_bc_vt = [1 ; 3];
    [F_opt,u_opt] = truss_solve(A,E,Truss_L,theta,elements,nof_nodes,F_lc,F_lc_hz,F_lc_vt, u_bc,u_bc_hz,u_bc_vt);
    node2disp(j) = u_opt(4);
end
[~,order] = min(abs(node2disp));
disp('4. Order of elements for minimum Node2 Veritcal Displacement')
truss_order = elements_switch(order,:) % ans = [1.5 2 2 1.5];
sprintf('4. Arrangements leads to minimum vertical displacement at node 2 are: \n 3 6 4 5 \n 6 3 4 5 \n 3 6 5 4 \n 6 3 5 4 \n From Left to Right')


%% User Defined Function
%For Concise, this function is used in Q4 to solve for force and displacements of this truss structure.
function [F,u] = truss_solve(A,E,Truss_L,theta,elements,nof_nodes,F_lc,F_lc_hz,F_lc_vt, u_bc,u_bc_hz,u_bc_vt)
k_coeff = A.*E./Truss_L;
nof_elements = size(elements,1);

% Prelocations
K = zeros(2*nof_nodes, 2*nof_nodes);
K_element = zeros(2*nof_nodes, 2*nof_nodes);
u = zeros(2*nof_nodes,1);
F = zeros(2*nof_nodes,1);

% Loading Condition
F_lc_hz = 2*F_lc_hz -1 ;
F_lc_vt = 2*F_lc_vt;
F_lc_index = [F_lc_hz ; F_lc_vt];
F_lc_index = sort(F_lc_index,'ascend');
F(F_lc_index) = F_lc;

% Boundary Conditions
u_bc_hz = 2*u_bc_hz -1 ;
u_bc_vt = 2*u_bc_vt ;
u_bc_index = [u_bc_hz ; u_bc_vt];
u_bc_index = sort(u_bc_index,'ascend');
u(u_bc_index) = u_bc;

% Construciton of K (Global)
for i = 1 : nof_elements
    R = [cos(theta(i)).^2 cos(theta(i))*sin(theta(i)); cos(theta(i))*sin(theta(i)) sin(theta(i)).^2];
    K_element(2*elements(i,1) - 1 : 2*elements(i,1) , 2*elements(i,1) - 1 : 2*elements(i,1)) = k_coeff(i)*R;
    K_element(2*elements(i,1) - 1 : 2*elements(i,1) , 2*elements(i,2) - 1 : 2*elements(i,2)) = k_coeff(i)*(-R);
    K_element(2*elements(i,2) - 1 : 2*elements(i,2) , 2*elements(i,1) - 1 : 2*elements(i,1)) = k_coeff(i)*(-R);
    K_element(2*elements(i,2) - 1 : 2*elements(i,2) , 2*elements(i,2) - 1 : 2*elements(i,2)) = k_coeff(i)*R;
    K = K + K_element;
    K_element(:,:) = 0;
end
% Nodal Displacements
node_index_all = transpose(1 : 2*nof_nodes);
F_known_index = setdiff(node_index_all,u_bc_index,'stable');
u(F_known_index) = K(F_known_index,F_known_index) \ F(F_known_index);

% Reactions Forces
F = K*u;

end