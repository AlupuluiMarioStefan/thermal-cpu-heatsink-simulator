close all;
clear all;

% Dimensions of the plate
xmin = 0; xmax = 2; ymin = 0; ymax = 3;

% Grid
Nx = 40; Ny = 60;
dltx = (xmax - xmin) / (Nx - 1);
dlty = (ymax - ymin) / (Ny - 1);
[x, y] = ndgrid(xmin + (0:Nx-1) * dltx, ymin + (0:Ny-1) * dlty);

% Source term
s = 15 * (x > 0.25 & x < 1.25 & y > 0.75 & y < 1.25);
figure; 
pcolor(x, y, s); 
xlabel('x'); ylabel('y'); title('Source');
colorbar; shading interp; colormap hot; axis equal; axis tight

% Build linear system for internal points
A = zeros(Nx * Ny, Nx * Ny);
b = zeros(Nx * Ny, 1);

% Internal grid points
for i = 2:Nx-1
    for j = 2:Ny-1
        C = (j - 1) * Nx + i; % Label of grid point (i,j)
        E = C + 1; % East neighbor
        W = C - 1; % West neighbor
        N = C + Nx; % North neighbor
        S = C - Nx; % South neighbor

        % Build the matrix (Laplacian)
        A(C, C) = -2 / dltx^2 - 2 / dlty^2;
        A(C, E) = 1 / dltx^2;
        A(C, W) = 1 / dltx^2;
        A(C, N) = 1 / dlty^2;
        A(C, S) = 1 / dlty^2;

        % Right-hand side
        b(C) = s(i, j);
    end
end

% Time stepping parameters
dltt = 0.1; % Time step
Nt = 30;    % Number of time steps

% Construct the implicit time-stepping matrix M = I - dltt * A
M = eye(Nx * Ny) - dltt * A;

% Apply Dirichlet Boundary Conditions (u = 0 on all boundaries)
% Top and Bottom boundaries
for i = 1:Nx
    % Bottom (j = 1)
    C_bot = i;
    M(C_bot, :) = 0; M(C_bot, C_bot) = 1;
    b(C_bot) = 0;
    
    % Top (j = Ny)
    C_top = (Ny - 1) * Nx + i;
    M(C_top, :) = 0; M(C_top, C_top) = 1;
    b(C_top) = 0;
end

% Left and Right boundaries
for j = 2:Ny-1
    % Left (i = 1)
    C_left = (j - 1) * Nx + 1;
    M(C_left, :) = 0; M(C_left, C_left) = 1;
    b(C_left) = 0;
    
    % Right (i = Nx)
    C_right = (j - 1) * Nx + Nx;
    M(C_right, :) = 0; M(C_right, C_right) = 1;
    b(C_right) = 0;
end

% Initial solution
u = zeros(Nx * Ny, 1); 

f = figure(2); 
hold on

for i = 1:Nt
    t = i * dltt;
    
    % Right-hand side for the implicit step: u^n + dltt * b
    rhs = u + dltt * b;
    
    % Explicitly zero out boundary entries in RHS to enforce u = 0
    for idx = 1:Nx
        rhs(idx) = 0;                    % Bottom boundary
        rhs((Ny - 1) * Nx + idx) = 0;    % Top boundary
    end
    for j = 2:Ny-1
        rhs((j - 1) * Nx + 1) = 0;       % Left boundary
        rhs((j - 1) * Nx + Nx) = 0;      % Right boundary
    end
    
    % Solve linear system for the next time step
    u = M \ rhs;
    
    % Plot solution
    um = reshape(u, Nx, Ny); 
    figure(2);
    pcolor(x, y, um); 
    shading interp; 
    colormap hot; 
    colorbar;
    xlabel('x'); ylabel('y'); 
    title(sprintf('t = %g', t));
    clim([0 3]); 
    axis equal; axis tight;
    
    M_movie(i) = getframe(f); 
end

figure(2);
movie(f, M_movie, 5);