close all; clear; clc

AE = 1/5;
L = pi/2; % Length of the whole domain
n = [5, 10, 20]; % Number of elements
u = zeros(21, 3); % u_1 is used to store the solution given by method 1
X = zeros(21,3);


for i = 1 : 3
    x = linspace(0, pi/2, n(i)+1).'; % Discretization
    X(1:n(i)+1,i) = x;
    l = L/n(i);
    k = AE/l; % "Spring Constant" in each local elements
    
    % Construction of K
    K = sparse(1:1+n(i), 1:1+n(i), [k, repmat(2*k,1, n(i)-1), k]) + ...
        sparse(1:n(i), 2:1+n(i), -k, n(i)+1, n(i)+1) + ...
        sparse(2:n(i)+1, 1:n(i), -k, n(i)+1, n(i)+1);
    
    % Construction of f
    x1 = x(1:end-1);
    x2 = x(2:end);
    F = l*[((x2 - x1).*cos(x1) + sin(x1) - sin(x2))./ ((x2 - x1).^2);0] + ...
        l*[0;((x1 - x2).*cos(x2) + sin(x2) - sin(x1))./ ((x2 - x1).^2)];
    
    % Solve Ku=F with BC W(0)=0
    u(2:n(i)+1,i) = K(2:n(i)+1, 2:n(i)+1)\F(2:n(i)+1); % Row 1 and Column1 eliminated by the
    
end

u

% Plots of FEM Solution
figure; hold on ; grid on ;
linestyle = ["bo-","k-.","r--"];
linewidth = [1,1,2];
for j = 1 : 3
    plot(X(1:n(j)+1,j),u(1:n(j)+1,j),linestyle(j),'LineWidth',linewidth(j))
end
title('FEM Solution With Various Elements'); xlabel('x'); ylabel('u')
legend('5 elements','10 elements','20 elements','Location','southeast'); hold off;

% Plots of Exaxt Solution
figure; grid on ;
x_exact = linspace(0,pi/2,20);
u_exact = 5*sin(x_exact);
plot(x_exact,u_exact,'g');
xlabel('x'); ylabel('u'); title('Exact Solution');

% Comparsion Between FEM and Exact Solution
figure; hold on; grid on;
for j = 1 : 3
    plot(X(1:n(j)+1,j),u(1:n(j)+1,j),linestyle(j),'LineWidth',linewidth(j))
end
plot(x_exact,u_exact,'g'); xlabel('x'); ylabel('u');
title('Comparsion Between FEM Solutions and Exact Solution');
legend('5 elements','10 elements','20 elements','Exact Solution','Location','southeast');

% Nov 19 2019 Late Turn in Part
u_x = zeros(20,3);
linestyle = ["bo-","k*-","rx-"];
figure; hold on;
for i = 1 : 3
    element_l = L/n(i); % Element Length
    u_x(1:n(i),i) = (1/element_l)*(u(2:n(i)+1,i) - u(1:n(i),i));
end

for i = 1 : 3
    x = linspace(0, pi/2, n(i)+1).';
    x_plot = repelem(x,2);
    x_plot(1) = []; x_plot(end) = [];
    u_x_plot = repelem(u_x(1:n(i),i),2);
    plot(x_plot,u_x_plot,linestyle(i));
%     for j = 1 : n(i)
%         x = linspace(0, pi/2, n(i)+1).';
%         plot(x(j:j+1), repmat(u_x(j,i),1,2), linestyle(i))
%     end
end

% Plots of Exaxt u_x
x_exact = linspace(0,pi/2,21);
u_x_exact = 5*cos(x_exact);
plot(x_exact,u_x_exact,'g');

xlabel('x')
ylabel('u,x','Rotation', 360)
title('FEM solution of u,x and Exact Solution');
legend('5 elements','10 elements','20 elements','Exact Solution')






