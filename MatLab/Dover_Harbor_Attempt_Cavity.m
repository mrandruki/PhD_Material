clc
clear

% Global parameters
global BOUNDARY
global FLOW
global SOLID
global FILENAME

BOUNDARY = 0;
FLOW = 1;
SOLID = 2;
FILENAME = 'Dover_Mean_Tidal_Flipped_OldBW';

% Define the domain
lenx = 1;                           % Domain length in the x direction
leny = 1;                           % Domain length in the y direction
nodes_per_unit = 400;
Nx = (lenx * nodes_per_unit);     % No. cells in the x direction
Ny = (leny * nodes_per_unit);     % No. cells in the y direction
dx = 1 / nodes_per_unit;            % x spacing
dy = 1 / nodes_per_unit;            % y spacing

% Model parameters
u0 = 0.1;                           % Velocity at top boundary
rho0 = 1;                           % Initial density
dt = 0.1;                          % Time step (what is the stability criteria?)
e = dx / dt;                        % Lattice Velocity
Re = 5000;                          % Reynolds number
tau = (6 * ((u0 * leny) / (e * dx * Re)) + 1) / 2;
nu = e * dx * (2 * tau - 1)/6;

% Solver parameters
iter = 0;                           % Iteration counter
max_iter = 15000;                   % Max no. iterations permitted
output_iter = 20;                   % No. iterations between solution output
time = 0;                           % Time variable
err = 1;                            % L^2 error between successive iterations
tol = 1d-4; %1d-5;                         % Convergence tolerance

w(9) = 4/9;
w(1:4) = 1/9;
w(5:8) = 1/36;
ex = [1 0 -1 0 1 -1 -1 1 0] * e;
ey = [0 1 0 -1 1 1 -1 -1 0] * e;
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];  % Opposite directions - not used

% Discretise the domain
x = (0 : Nx) * dx;
y = (0 : Ny) * dy;

% Initalise the flow, solid and boundaries cells
% material_d = FLOW * ones(Nx, Ny);

% Load the material from a file
load(strcat(FILENAME, '.mat'));
% material(360:375,:) = [];
material_d = material;
% material_d(1 : Nx, 240 : Ny) = 1;

% Initialise solution arrays
f = zeros(Nx, Ny, 9);
feq = zeros(Nx, Ny, 9);
ftemp = zeros(Nx, Ny, 9);
rho = rho0 * ones(Nx, Ny);
u = zeros(Nx, Ny);
v = zeros(Nx, Ny);
rho_old = rho;
u(material_d(Nx, 1) == FLOW) = u0;

% Main Loop
tic
while err > tol % Perform iterations until convergence % for iter = 1 : max_iter
    
    % Compute Feq
    t1 = (u.*u + v.*v) / (e*e);
    for a = 1 : 9
        t2 = (u * ex(a) + v * ey(a)) / (e*e);
        feq(:, :, a) = rho * w(a) .* (1 + 3*t2 + 4.5*t2.*t2 - 1.5*t1);
    end
    
    % Collide Stream
    for y = 1 : Ny
        yp = y + 1;
        yn = y - 1;
        for x = 1 : Nx
            xp = x + 1;
            xn = x - 1;
            if (material_d(x,y) == BOUNDARY || material_d(x,y) == FLOW)
                if xp <= Nx && (material_d(xp,y) == BOUNDARY || material_d(xp,y) == FLOW)
                    ftemp(xp,y,1) = f(x,y,1) - (f(x,y,1) - feq(x,y,1)) / tau;
                end
                if yp <= Ny && (material_d(x,yp) == BOUNDARY || material_d(x,yp) == FLOW)
                    ftemp(x,yp,2) = f(x,y,2) - (f(x,y,2) - feq(x,y,2)) / tau;
                end
                if xn >= 1 && (material_d(xn,y) == BOUNDARY || material_d(xn,y) == FLOW)
                    ftemp(xn,y,3) = f(x,y,3) - (f(x,y,3) - feq(x,y,3)) / tau;
                end
                if yn >= 1 && (material_d(x,yn) == BOUNDARY || material_d(x,yn) == FLOW)
                    ftemp(x,yn,4) = f(x,y,4) - (f(x,y,4) - feq(x,y,4)) / tau;
                end
                if xp <= Nx && yp <= Ny && (material_d(xp,yp) == BOUNDARY || material_d(xp,yp) == FLOW)
                    ftemp(xp,yp,5) = f(x,y,5) - (f(x,y,5) - feq(x,y,5)) / tau;
                end
                if xn >= 1 && yp <= Ny && (material_d(xn,yp) == BOUNDARY || material_d(xn,yp) == FLOW)
                    ftemp(xn,yp,6) = f(x,y,6) - (f(x,y,6) - feq(x,y,6)) / tau;
                end
                if xn >= 1 && yn >= 1 && (material_d(xn,yn) == BOUNDARY || material_d(xn,yn) == FLOW)
                    ftemp(xn,yn,7) = f(x,y,7) - (f(x,y,7) - feq(x,y,7)) / tau;
                end
                if xp <= Nx && yn >= 1 && (material_d(xp,yn) == BOUNDARY || material_d(xp,yn) == FLOW)
                    ftemp(xp,yn,8) = f(x,y,8) - (f(x,y,8) - feq(x,y,8)) / tau;
                end
                ftemp(x,y,9) = f(x,y,9) - (f(x,y,9) - feq(x,y,9)) / tau;
            else
                ftemp(x,y,:) = 0;
            end
        end
    end
    
    % Bounceback
    for y = 2 : Ny - 1
        yp = y + 1;
        yn = y - 1;
        for x = 2 : Nx - 1
            xp = x + 1;
            xn = x - 1;
            if material_d(x,y) == FLOW
                if material_d(xp,y) == SOLID
                    ftemp(x,y,3) = ftemp(x,y,1);
                end
                if material_d(x,yp) == SOLID
                    ftemp(x,y,4) = ftemp(x,y,2);
                end
                if material_d(xn,y) == SOLID
                    ftemp(x,y,1) = ftemp(x,y,3);
                end
                if material_d(x,yn) == SOLID
                    ftemp(x,y,2) = ftemp(x,y,4);
                end
                if material_d(xp,yp) == SOLID
                    ftemp(x,y,7) = ftemp(x,y,5);
                end
                if material_d(xn,yp) == SOLID
                    ftemp(x,y,8) = ftemp(x,y,6);
                end
                if material_d(xn,yn) == SOLID
                    ftemp(x,y,5) = ftemp(x,y,7);
                end
                if material_d(xp,yn) == SOLID
                    ftemp(x,y,6) = ftemp(x,y,8);
                end
            end
        end
    end
    
    % Boundaries
    % Outlet along the east boundary
%     ftemp(Nx, :, :) = ftemp(Nx - 1, :, :);
%     u(Nx, :) = u(Nx - 1, :);
%     v(Nx, :) = v(Nx - 1, :);
%     rho(Nx, :) = rho(Nx - 1, :);
% 
%     % Bounce along the east boundary
    ftemp(Nx, :, 3) = ftemp(Nx, :, 1);
    ftemp(Nx, :, 7) = ftemp(Nx, :, 5);
    ftemp(Nx, :, 6) = ftemp(Nx, :, 8);
    
    % Bounce along the west boundary
    ftemp(1, :, 1) = ftemp(1, :, 3);
    ftemp(1, :, 5) = ftemp(1, :, 7);
    ftemp(1, :, 8) = ftemp(1, :, 6);
    
    % Bounce along the south boundary
    ftemp(:, 1, 2) = ftemp(:, 1, 4);
    ftemp(:, 1, 5) = ftemp(:, 1, 7);
    ftemp(:, 1, 6) = ftemp(:, 1, 8);
    
    % Across the lid
    rhon = ftemp(2:Nx-1, Ny, 9) + ftemp(2:Nx-1, Ny, 1) + ftemp(2:Nx-1, Ny, 3) + ...
         2 * (ftemp(2:Nx-1, Ny, 2) + ftemp(2:Nx-1, Ny, 5) + ftemp(2:Nx-1, Ny, 6));
    ftemp(2:Nx-1, Ny, 4) = ftemp(2:Nx-1, Ny, 1) + ftemp(2:Nx-1, Ny, 3) + ftemp(2:Nx-1, Ny, 2) + ...
        2 * (ftemp(2:Nx-1, Ny, 5) + ftemp(2:Nx-1, Ny, 6)) - rhon/3 - rhon * u0^2;
    ftemp(2:Nx-1, Ny, 7) = rhon/6 + 1/2 * rhon * u0^2 - 1/2 * rhon * u0 - ...
        ftemp(2:Nx-1, Ny, 3) - ftemp(2:Nx-1, Ny, 6);
    ftemp(2:Nx-1, Ny, 8) = rhon/6 + 1/2 * rhon * u0^2 + 1/2 * rhon * u0 - ...
        ftemp(2:Nx-1, Ny, 1) - ftemp(2:Nx-1, Ny, 5);
    
    % Northwest corner
    rhonw = ftemp(1, Ny, 9) + 2 * ftemp(1, Ny, 3) + 4 * ftemp(1, Ny, 6) + ...
         2 * ftemp(1, Ny, 2);
     ftemp(1, Ny, 1) = 2 * rhonw / 3 - ftemp(1, Ny, 9) - ftemp(1, Ny, 3);
     ftemp(1, Ny, 4) = 2 * rhonw / 3 - ftemp(1, Ny, 9) - ftemp(1, Ny, 2);
     ftemp(1, Ny, 5) = rhonw / 6 - ftemp(1, Ny, 2) - ftemp(1, Ny, 6);
     ftemp(1, Ny, 7) = rhonw / 6 - ftemp(1, Ny, 3) - ftemp(1, Ny, 6);
     ftemp(1, Ny, 8) = -2 * rhonw / 3 + ftemp(1, Ny, 9) + ftemp(1, Ny, 2) + ...
         ftemp(1, Ny, 3) + f(1, Ny, 6);
     
     % Northeast corner
    rhone = ftemp(1, Ny, 9) + 2 * ftemp(1, Ny, 1) + 4 * ftemp(1, Ny, 5) + ...
         2 * ftemp(1, Ny, 2);
     ftemp(1, Ny, 3) = 2 * rhone / 3 - ftemp(1, Ny, 9) - ftemp(1, Ny, 1);
     ftemp(1, Ny, 4) = 2 * rhone / 3 - ftemp(1, Ny, 9) - ftemp(1, Ny, 2);
     ftemp(1, Ny, 6) = rhone / 6 - ftemp(1, Ny, 2) - ftemp(1, Ny, 5);
     ftemp(1, Ny, 8) = rhone / 6 - ftemp(1, Ny, 1) - ftemp(1, Ny, 5);
     ftemp(1, Ny, 7) = -2 * rhone / 3 + ftemp(1, Ny, 9) + ftemp(1, Ny, 2) + ...
         ftemp(1, Ny, 3) + ftemp(1, Ny, 6);

    % Calculating the components
    f = ftemp;
    sum_u = zeros(size(u));
    sum_v = zeros(size(v));
    sum_rho = zeros(size(rho));
    flowcells = material_d == FLOW;
    for a = 1 : 9
        sum_rho = sum_rho + f(:, :, a) .* flowcells;
        sum_u = sum_u + ex(a) * f(:, :, a) .* flowcells;
        sum_v = sum_v + ey(a) * f(:, :, a) .* flowcells;
    end
    rho = sum_rho;
    u = sum_u ./ rho;
    v = sum_v ./ rho;
    
    % Clear solid cells
    u(material_d == SOLID) = 0;
    v(material_d == SOLID) = 0;
    
    % Update progress indicators
    iter = iter + 1;
    time = iter * dt;
    err = sqrt(sum(sum((rho - rho_old).^2)));
    rho_old = rho;
    
    % Exit loop it max iteration number is exceeded
    if iter == max_iter
        fprintf('Warning! Max iterations exceeded.\n');
        cont = input('Continue with the calculations? (y/n) > ', 's');
        save(strcat(FILENAME, '_Outflow_Re', int2str(Re)));
        if cont == 'y'
            max_iter = max_iter + 25000;
        else
            break
        end
    end
    
%     Exit loop if error is large
    if err > 10000
        fprintf('Warning! Error = %1.2e. Terminating program.\n', err);
        break
    end
    
    % Output solution
    if mod(iter, output_iter) == 0
        
        % Print to command window
        fprintf('iteration : %6i | error : %1.2e\n', iter, err);
        
        % Plot solution
        nx = 5;
        ny = 5;
        xr = (0:1:Nx-1)*dx;
        yr = (0:1:Ny-1)*dy;
        [xc,yc] = ndgrid(xr,yr);
        quiver(xc(1:nx:Nx, 1:ny:Ny), yc(1:nx:Nx, 1:ny:Ny), u(1:nx:Nx, 1:ny:Ny), v(1:nx:Nx, 1:ny:Ny),3);
        axis([0 lenx 0 leny]);
        drawnow  
    end
end
cpu_time = toc;
save(strcat(FILENAME, '_Outflow_Re', int2str(Re)));

% Output final information to command window
fprintf('-------------------------------------\n')
fprintf('iterations : %6i\n', iter)
fprintf('CPU time   : %8.3fs\n', cpu_time)

% ------------------------------------------------------------------------
function [i, j, ii, jj] = get_cell_indices(dir, Nx, Ny)

% This function returns the indices of the cells which are under 
% consideration (i and j) and the neighbouring cell determined by the 
% direction dir (ii and jj).

if dir == 1 % East
    i = 1 : Nx - 1;
    j = 1 : Ny;
    ii = i + 1;
    jj = j;
elseif dir == 2 % North
    i = 1 : Nx;
    j = 1 : Ny - 1;
    ii = i;
    jj = j + 1;
elseif dir == 3 % West
    i = 2 : Nx;
    j = 1 : Ny;
    ii = i - 1;
    jj = j;
elseif dir == 4 % South
    i = 1 : Nx;
    j = 2 : Ny;
    ii = i;
    jj = 1 : Ny - 1;
elseif dir == 5 % Northeast
    i = 1 : Nx - 1;
    j = 1 : Ny - 1;
    ii = i + 1;
    jj = j + 1;
elseif dir == 6 % Northwest
    i = 2 : Nx;
    j = 1 : Ny - 1;
    ii = 1 : Nx - 1;
    jj = 2 : Ny;
elseif dir == 7 % Southwest
    i = 2 : Nx;
    j = 2 : Ny;
    ii = i - 1;
    jj = j - 1;
elseif dir == 8 % Southeast
    i = 1 : Nx-1;
    j = 2 : Ny;
    ii = i + 1;
    jj = j - 1;
elseif dir == 9 % Centre
    i = 1 : Nx;
    j = 1 : Ny;
    ii = i;
    jj = j;
end

end
