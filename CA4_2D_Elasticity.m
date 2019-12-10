close all;clear;clc
% Physical Parameters
k = 50;

nof_elem = [2,4,8];

% Shape Function
N = @(xi, eta) 0.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1-xi)*(1+eta) (1+xi)*(1+eta)];

% Integrand of Element Stiffness Matrix
B = @(xi,eta,inv_A) inv_A*0.25*[eta-1, 1-eta, -1-eta, 1+eta; xi-1, -1-xi, 1-xi, 1+xi];

% Body Force
Q = @(vec_x, vec_y, N) sin(0.5*pi*N*vec_x)*sin(0.5*pi*N*vec_y); % Isoparametric Mapping; x and y are position vector of the element

% Gauss Integration Points
xi_1 = -1/sqrt(3); xi_2 = 1/sqrt(3);
eta_1 = -1/sqrt(3); eta_2 = 1/sqrt(3);
gip = [xi_1 eta_1; xi_1 eta_2; xi_2 eta_1; xi_2 eta_2];

% % Jacobian
% Jv = @(xi,eta,vec_x,vec_y)

linestyle = ["bo-","k-.","r--"];

u_all = cell(3,1);
x_all = cell(3,1);
u_x_all = cell(3,1);
error_all = cell(3,1);
x_exact_all = cell(3,1);
u_c_all = zeros(3,1);
h = zeros(3,1);

for nelem = 1 : 3
    
    % Domain Discretization
    nx = nof_elem(nelem); ny = nof_elem(nelem); % Number of elements in x and y directions
    Lx = 1 ; Ly = 1; % Domain length in x and y directions
    dx = Lx/nx; dy = Ly/ny; % element length for each elements in x and y directions
    h(nelem) = dx;
    grid = reshape(linspace(1,(nx+1)*(ny+1),(nx+1)*(ny+1)),ny+1,nx+1).'; % Nodal Number
    
    % Jacobians
    Jv = 0.25*dx*dy;
    
    % A Inverse
    inv_A = [2/dx 0; 0 2/dy];
    
    % Nodal coordinate values for all elements; Stored by Rows; size: (nx*ny,4)
    X = repmat([linspace(0,Lx-Lx/nx,nx).' linspace(Lx/nx,Lx,nx).'],ny,2);
    Y = [repmat(repelem(linspace(0,Ly-Ly/ny,ny),1,nx).',1,2) repmat(repelem(linspace(Ly/ny,Ly,ny),1,nx).',1,2)];
    
    %% Mapping Table
    a = 1:1:(nx+1)*ny;
    a(mod(a,nx+1)==0) = [];
    b = setdiff(1:1:(nx+1)*ny,ones(1,ny)+(nx+1)*linspace(0,ny-1,ny));
    index = ones(4*nx*ny,1);
    index(1:4:4*nx*ny,1) = a;
    index(2:4:4*nx*ny,1) = b;
    index(3:4:4*nx*ny,1) = a + (nx+1);
    index(4:4:4*nx*ny,1) = b + (nx+1);
    index = reshape(index,4,size(index,1)/4).';
    
    %% Construction of Stiffness Matrix
    K = zeros((nx+1)*(ny+1),(nx+1)*(ny+1));
    f = zeros((nx+1)*(ny+1),1);
    
    for i = 1 : nx*ny % Scanning through all elements
        
        K_e = zeros(4,4);
        f_e = zeros(4,1);
        
        % Gauss Integration of BTB
        for j = 1 : 4
            N_e = N(gip(j,1),gip(j,2));
            Q_e = Q(X(i,:).', Y(i,:).', N_e);
            B_e = B(gip(j,1), gip(j,2), inv_A);
            
            K_e = K_e + B_e.'*B_e;
            f_e = f_e + N_e.'*Q_e;
        end
        
        K(index(i,:),index(i,:)) = K(index(i,:),index(i,:)) + k*Jv*K_e;
        f(index(i,:)) = f(index(i,:)) + Jv*f_e;
        
    end
    
    % Apply Boundary Condition and solve the system of equations
    % u = 0 at boundary Gamma1 and Gamma2
    node_re = setdiff(grid, [grid(1,:), grid(:,1).']);
    grid_re = reshape(node_re,nx,ny).';
    u_re = K(node_re,node_re)\f(node_re);
    
    % Q1. u along line y = 0.5
    u_re = reshape(u_re,nx,ny).';
    u = [0 u_re(ny/2, :)]; % u along the line of y = 0.5
    u_all{nelem} = u;
    x = linspace(0,Lx,nx+1); x_all{nelem} = x;
    
    % Q2. dudx along line y = 0.5
    G = [-1/dx, 1/dx];
    u_x = G*[u(1:end-1); u(2:end)]; u_x_all{nelem} = u_x;

    % Q3. error along line y = 0.5 (Interpolation Applied)
    x_exact = linspace(0,1,10*nof_elem(nelem));
    u_exact = 2/(50*pi^2)*sin(0.5*pi*x_exact).*sin(0.5*pi*0.5);
    u_error = zeros(10*nof_elem(nelem),1);
    xi = linspace(-1,1,10); % 10 interpolation points for each element
    for j = 1 : nof_elem(nelem)
        d1 = u(1:end-1); d2 = u(2:end);
        u_interp = 0.5*(d1(j) + d2(j)) + 0.5*(d2(j) - d1(j))*xi;
        u_error(1+(j-1)*10:j*10) = abs(u_exact(1+(j-1)*10:10*j) - u_interp);
    end    
    error_all{nelem} = u_error;
    x_exact_all{nelem} = x_exact;

    % Q4. error at point (0.5,0.5); u_c = u @ Domain Center
    u_c = u(nx/2+1) - 2/(50*pi^2)*sin(0.5*pi*0.5)*sin(0.5*pi*0.5); u_c_all(nelem) = u_c;
    
end

%% Q1. Plot of u
figure; hold on ;
for i = 1 : 3
    plot(x_all{i}, u_all{i},linestyle(i));
end
x_exact = linspace(0,1,100);
u_exact = 2/(50*pi^2)*sin(0.5*pi*x_exact).*sin(0.5*pi*0.5); % y = 0.5
plot(x_exact, u_exact,'g');
hold off

lgd1 = legend('$$2\times2 FEM$$','$$4\times4 FEM$$','$$8\times8 FEM$$','Exact Solution','Location','southeast');
set(lgd1,'Interpreter','latex');
xlabel('x');
ylabel('u_{(x,0.5)}');
title('FEM Solution and Exact Solution of u^h')

%% Q2. Plot of u_x
figure; hold on;
for i = 1 : 3
    x_plot = x_all{i};
    x_plot = repelem(x_plot,2);
    x_plot(1) = []; x_plot(end) = [];
    u_x_plot = u_x_all{i};
    u_x_plot = repelem(u_x_plot(1:nof_elem(i)),2);
    plot(x_plot,u_x_plot,linestyle(i));
end
dudx_exact = 1/(50*pi)*cos(0.5*pi*x_exact).*sin(0.5*pi*0.5);
plot(x_exact, dudx_exact)
hold off; 
xlabel('x');
ylabel('u^h_{,x}')
title('FEM and Exact Solution of u^h_{,x}')
lgd2 = legend('$$2\times2 FEM$$','$$4\times4 FEM$$','$$8\times8 FEM$$','Exact Solution','Location','northeast');
set(lgd2,'Interpreter','latex');

%% Q3. Plot of Error
figure; hold on ;
for i = 1 : 3
    plot(x_exact_all{i}, error_all{i},linestyle(i));
end
hold off; 
xlabel('x');
yl = ylabel('$$|u-u^{h}|$$');
set(yl,'Interpreter','Latex');
lgd3 = legend('$$2\times2 FEM$$','$$4\times4 FEM$$','$$8\times8 FEM$$','Location','northeast');
set(lgd3,'Interpreter','latex');
title('Error of FEM Solution')

%% Plot of Convergence Rate
figure
loglog(h, u_c_all,'-o');
xlabel('Element Size h');
yl = ylabel('$$|u-u^{h}|$$');
set(yl,'Interpreter','Latex');
title('Rate of Convergence of Mesh Refinement')



