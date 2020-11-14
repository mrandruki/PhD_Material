clear 
clc

global Filename
Filename = 'Dover_Mean_Tidal_Flipped.mat';

length = 20;
width = 16;
h0 = 1.5;
v0 = 0;
Q0 = 2.2;
Lx = 201;
Ly = 161;
max_t = 1;            
gacl = 9.81;
fb = 0.00248; % fb=g/((h^1/6)/nb)^2, nb=0.012% Manning's coefficient
D=0.5;

% Shape of the support beams (shape type, x, y, radius)
 shapes = [];

% Changes the shape of the bed (y, h)
bed_shape = [1 0.333; 20 0.2; 40 0.24; 80 0.28; 120 0.32; 161 0.35;];

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
trapezium1 =[5 1.6; 6 1.4; 8 1.6; 7 1.4;];
trapezium2 =[6 1.8; 5 1.6; 7 1.8; 8 1.6;];
trapezium3 =[5 3.2; 6 3; 8 3.2; 7 3;];
trapezium4 =[6 3.4; 5 3.2; 7 3.4; 8 3.2;];
trapezium5 =[5 4.8; 6 4.8; 8 4.8; 7 4.6;];
trapezium6 =[6 5; 5 4.8; 7 5; 8 4.8;];
trapezium7 =[5 6.4; 6 6.2; 8 6.4; 7 6.2;];
trapezium8 =[6 6.6; 5 6.4; 7 6.6; 8 6.4;];
trapezium9 =[5 8; 6 7.8; 8 8; 7 7.8;];
trapezium10 =[6 8.2; 5 8; 7 8.2; 8 8;];
trapezium11 =[5 9.6; 6 9.4; 8 9.6; 7 9.4;];
trapezium12 =[6 9.8; 5 9.6; 7 9.8; 8 9.6;];
trapezium13 =[5 11.2; 6 11; 8 11.2; 7 11;];
trapezium14 =[6 11.4; 5 11.2; 7 11.4; 8 11.2;];
trapezium15 =[5 12.8; 6 12.6; 8 12.8; 7 12.6;];
trapezium16 =[6 13; 5 12.8; 7 13; 8 12.8;];
trapezium17 =[5 14.4; 6 14.2; 8 14.4; 7 14.2;];
trapezium18 =[6 14.6; 5 14.4; 7 14.6; 8 14.4;];

% Determine solid values
solid = solid_values(Lx,Ly,dx,dy,shapes);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium1);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium2);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium3);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium4);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium5);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium6);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium7);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium8);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium9);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium10);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium11);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium12);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium13);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium14);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium15);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium16);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium17);
solid = trapezium_values(Lx,Ly,dx,dy,solid,trapezium18);

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


% Main loop
iteration=0;
time=0; 
    
while iteration < max_t

    iteration=iteration+1;
    
    time=iteration*dt;
        
    feq = compute_feq(Lx,Ly,h,u,v,e,ex,ey,gacl,solid);          
    
    ftemp = collide_stream(Xs,Xe,Ys,Ye,tau,gacl,dt,Sx,Sy,u,v,e,ex,ey,feq,f,h,fb,solid); 

    ftemp=BC_Body_slipNS(Xs,Xe,Ly,ftemp);
    
    ftemp=BC_Body_main(Xs,Xe,Ys,Ye,solid,ftemp);
    
    ftemp=BC_InOut(Lx,ftemp);
    
    [h,u,v,f]=solution(Xs,Xe,Ys,Ye,solid,ftemp,ex,ey,f); 

    [h,u,v]=BC_InflowOutflow(h,u,v,Q0,dy,h0,Lx,Ly); 
     
    disp(sprintf('%d  %f', iteration, h(95,51)));
   
    if iteration == max_t

        %write file
        Mat=write_file(Lx,Ly,dx,dy,dt,u,v,h,zb,tau, iteration, length, width, h0, max_t, nu, fb, Re,Fr);
        save('result');
        
        count = input('Contiue computations (y/n)?', 's');
        
        if count == 'y';
            addItea = input('Type additional number: ');
            max_t = max_t+addItea;            
        end
    end
end

[Q]=continuity(Xs,Xe,Ys,Ye,h,u,dy);   