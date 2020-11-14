% Plot streamlines

% The numerical data
load Dover_Mean_Tidal_Flipped_NewBW_SWE_LW_iter_500.mat

[x,y]=meshgrid(1:Lx,1:Ly);
streamslice((x-1)/(Lx-1),(y-1)/(Ly-1),squeeze(u)',squeeze(v)',3);

% [x,y]=meshgrid(1:Nx,1:Ny);
% streamslice((x-1)/(Nx-1),(y-1)/(Ny-1),squeeze(u)',squeeze(v)',5);

%   u = reshape(sqrt(ux.^2+uy.^2),Lx,Ly);
%     u(bbRegion) = nan;
%     imagesc(u(:,Ly:-1:1)'./uLid);
%     colorbar
%     axis equal off; drawnow
% 
% figure,contour(stream',60,'k-')
% 
% % title ('Streamline of Lid-Driven Cavity at Re 1000 with a cuboid within the centre')
% xlabel('x axis')
% ylabel('y axis')

% Plot solution
%         nx = 10;
%         ny = 10;
%         xr = (0:1:Lx-1)*dx;
%         yr = (0:1:Ly-1)*dy;
%         [xc,yc] = ndgrid(xr,yr);
%         quiver(xc(1:nx:Lx, 1:ny:Ly), yc(1:nx:Lx, 1:ny:Ly), u(1:nx:Lx, 1:ny:Ly), v(1:nx:Lx, 1:ny:Ly),2);
%         axis([0 length 0 width]);
%         drawnow  
