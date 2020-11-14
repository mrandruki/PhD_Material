clear 
clc

global Filename
Filename = 'Dover_Mean_Tidal_Flipped_NewBW';

load(strcat(Filename, '.mat'));

length = 1;
width = 1;
h0 = 4;
v0 = 0;
Q0 = 2.2;
Lx = 400;
Ly = 400;         
gacl = 9.81;
fb = 0.00248; % fb=g/((h^1/6)/nb)^2, nb=0.012% Manning's coefficient
D=0.5;

% Shape of the support beams (shape type, x, y, radius)
 shapes = [];

% Changes the shape of the bed (y, h)
bed_shape = [1 0.1; 399 0.1;];

% Calculation of parameters
dx = length/(Lx-1);    % dx=length/(Lx-1)
u0 = 0.5;    % u0=Q0/(h0*width)
dy = dx;
e = 10.;
dt = dx/e;
tau = 0.5*(1.+1.E-02*6.*dt/(dx*dx));
nu = e*dx*(2*tau-1)/6;
Xs = 1;Xe = Lx;Ys = 1;Ye = Ly; 
nermax = (Lx-1)*(Ly-1);
Fr = u0/sqrt(gacl*h0);
Re = u0*h0/nu;
ReD = u0*D/nu;

% Shape of the bank supports (top left, bottom left, top right, bottom right)

% Determine solid values
% solid = solid_values(Lx,Ly,dx,dy,shapes);
solid = material;

% Bed data and slope calculation
zb = bed_data(Lx,Ly,dx,bed_shape);
[Sx,Sy]=bed_slope(zb,Lx,Ly,dx,dy);

% Initialize the depth and velocity field
for y=1:Ly
    for x=1:Lx
        if solid(x,y)<=1
            h(x,y) = h0-zb(x,y); 
            v(x,y) = 0;
            u(x,y) = u0;
        else
            h(x,y) = 0;
            u(x,y) = 0;
            v(x,y) = 0;
        end
    end
end

[ex,ey]=setup(e);                                   % calculate ex and ey
feq = compute_feq(Lx,Ly,h,u,v,e,ex,ey,gacl,solid);  % initial feq
f=feq;                                              % set initial feq to f


% Solver parameters
iteration = 0;                           % Iteration counter
max_t = 5000;                   % Max no. iterations permitted
output_iter = 50;                   % No. iterations between solution output
time = 0;                           % Time variable
err = 1;                            % L^2 error between successive iterations
tol = 1d-4; %1d-5;                         % Convergence tolerance
    
while iteration < max_t

    iteration=iteration+1;
    
    time=iteration*dt;
        
    feq = compute_feq(Lx,Ly,h,u,v,e,ex,ey,gacl,solid);          
    
    ftemp = collide_stream(Xs,Xe,Ys,Ye,tau,gacl,dt,Sx,Sy,u,v,e,ex,ey,feq,f,h,fb,solid); 

    ftemp=BC_Body_slipNS(Xs,Xe,Ly,ftemp);
    
    ftemp=BC_Body_main(Xs,Xe,Ys,Ye,solid,ftemp);
    
    ftemp=BC_InOut(Lx,Ly,ftemp);
    
    [h,u,v,f]=solution(Xs,Xe,Ys,Ye,solid,ftemp,ex,ey,f); 

    [h,u,v]=BC_InflowOutflow(h,u,v,Q0,dy,dx,h0,u0,Lx,Ly); 
    
    if mod(iteration, output_iter) == 0
        
        % Print to command window
        fprintf('%d  %f\n', iteration, h(round(Lx / 2), round(Ly / 2)));
        
        % Plot solution
        nx = 5;
        ny = 5;
        xr = (0:1:Lx-1)*dx;
        yr = (0:1:Ly-1)*dy;
        [xc,yc] = ndgrid(xr,yr);
        quiver(xc(1:nx:Lx, 1:ny:Ly), yc(1:nx:Lx, 1:ny:Ly), u(1:nx:Lx, 1:ny:Ly), v(1:nx:Lx, 1:ny:Ly),3);
        axis([0 length 0 width]);
        drawnow  
    end
    
    if mod(iteration, 500) == 0
        save(strcat(Filename, '_SWE_LW_iter_', int2str(iteration)));
    end
   
    if iteration == max_t

        %write file
%         Mat=write_file(Lx,Ly,dx,dy,dt,u,v,h,zb,tau, iteration, length, width, h0, max_t, nu, fb, Re,Fr);
%         save(strcat(Filename, '_SWE'));
        
        count = input('Contiue computations (y/n)?', 's');
        
        if count == 'y';
            addItea = input('Type additional number: ');
            max_t = max_t+addItea;            
        end
    end
end

[Q]=continuity(Xs,Xe,Ys,Ye,h,u,dy);   