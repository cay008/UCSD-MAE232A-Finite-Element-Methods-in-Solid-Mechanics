close all;clear;clc
% Physical Parameters
k = 50; Jv = 1/16;

% Domain Discretization
nx = 2; ny = 2; Lx = 1 ; Ly = 1;
grid = reshape(linspace(1,(nx+1)*(ny+1),(nx+1)*(ny+1)),ny+1,nx+1).';

% Find A inverse for all elements
x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1);
x = repmat(x.',1,nx); % [x1 x1 ... x1 ; x2 x2 ... x2 ; ...... ; x_ny x_ny ... x_ny]; Size : nx+1 by ny
y = repmat(y,ny,1);   % [y1 y2 ... y_ny ; y1 y2 ... y_ny ; ...... ; y1 y2 ... y_ny]; Size : nx+1 by ny+1

x = x(2:nx+1,:) - x(1:nx,:); % [x2-x1 x2-x1 ... x2-x1 ; x3-x2 x3-x2 ... x3-x2 ; ...... ; x_n+1-x_n x_n+1-x_n ... x_n+1-x_n]
y = y(:,2:ny+1) - y(:,1:ny); % [y2-y1 y3-y2 ... y_n-y_n-1; repeat; repeat; .....; repeat];

x = 2./reshape(x,[],1); 
x = [x;0];
y = 2./reshape(y,[],1); % size: nof_element by 1
y = [0;y]; 

inv_A = sparse([1:2:2*nx*ny 2*nx*ny],[1:2:2*nx*ny 2*nx*ny],x) + sparse([1 2:2:2*nx*ny], [1 2:2:2*nx*ny],y);

% Anonymous Functions
BTB = @(xi,eta,inv_A) transpose(inv_A*1/4*[eta-1 1-eta -1-eta 1+eta; xi-1 -1-xi 1-xi 1+xi])*...
    inv_A*1/4*[eta-1 1-eta -1-eta 1+eta; xi-1 -1-xi 1-xi 1+xi];

NTQ = @(xi,eta,x,y) sin(0.5*pi*1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1-xi)*(1+eta) (1+xi)*(1+eta)]*x)*...
    sin(0.5*pi*1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1-xi)*(1+eta) (1+xi)*(1+eta)]*y)*...
    1/4*[(1-xi)*(1-eta); (1+xi)*(1-eta); (1-xi)*(1+eta); (1+xi)*(1+eta)];

% Gauss Integration Points
xi_1 = -1/sqrt(3); xi_2 = 1/sqrt(3);
eta_1 = xi_1; eta_2 = xi_2;
gip = [xi_1 eta_1; xi_1 eta_2; xi_2 eta_1; xi_2 eta_2];

% x and y coordinates used for all elements
X = repmat([linspace(0,Lx-Lx/nx,nx).' linspace(Lx/nx,Lx,nx).'],ny,2);
Y = [repmat(repelem(linspace(0,Ly-Ly/ny,ny),1,nx).',1,2) repmat(repelem(linspace(Ly/ny,Ly,ny),1,nx).',1,2)];

% Prelocations
K = zeros(nx*ny,nx*ny);
f = zeros(nx*ny,1);

% Where to assemble the local stiffness matrices?
a = 1:1:(nx+1)*ny;
a(mod(a,nx+1)==0) = [];
b = setdiff(1:1:(nx+1)*ny,ones(1,ny)+(nx+1)*linspace(0,ny-1,ny));
index = ones(4*nx*ny,1);
index(1:4:4*nx*ny,1) = a;
index(2:4:4*nx*ny,1) = b;
index(3:4:4*nx*ny,1) = a + (nx+1);
index(4:4:4*nx*ny,1) = b + (nx+1);
index = reshape(index,4,size(index,1)/4).';

% Construction of Stiffness Matrix
for i = 1 : nx*ny % Scanning through all elements
    
    K_e = zeros(4,4);
    f_e = zeros(4,1);
    
    % Gauss Integration of BTB
    for j = 1 : 4
    K_e = K_e + BTB(gip(j,1),gip(j,2),inv_A(2*i-1:2*i,2*i-1:2*i));
    f_e = f_e + NTQ(gip(j,1),gip(j,2),X(i,:).',Y(i,:).');
    end
    
    % Display element stiffness matrix and force vector of element 1
    % (following the index given in CA3);
    if i == 4
        K_1 = K_e;
        K_1(:,3:4) = [K_1(:,4) K_1(:,3)];
        K_1(3:4,:) = [K_1(4,:); K_1(3,:)];
        k*Jv*K_1
        
        f_1 = f_e;
        f_1(3:4) = [f_1(4);f_1(3)];
        Jv*f_1
    end

    K(index(i,:),index(i,:)) = k*Jv*K_e;
    f(index(i,:)) = Jv*f_e;
    
end

    
    

