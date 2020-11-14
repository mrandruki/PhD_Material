% Plot streamlines

% The numerical data
load Dover_Mean_Tidal_Flipped_OldBW_Outflow_Re5000.mat

% stream = zeros(Nx,Ny);
%   for i = 1 : Nx 
%      stream(i,1) = 0.;
%      for j = 2 : Ny
%         stream(i,j) = stream(i,j-1) + 0.5*(rho(i,j)*u(i,j)+rho(i,j-1)*u(i,j-1))*dy;
%      end
%   end

% [x,y]=meshgrid(1:lx,1:ly);
% streamslice((x-1)/(lx-1),(y-1)/(ly-1),squeeze(ux)',squeeze(uy)',1);
% 
[x,y]=meshgrid(1:Nx,1:Ny);
streamslice((x-1)/(Nx-1),(y-1)/(Ny-1),squeeze(u)',squeeze(v)',1);

%   u = reshape(sqrt(ux.^2+uy.^2),lx,ly);
%     u(bbRegion) = nan;
%     imagesc(u(:,ly:-1:1)'./uLid);
%     colorbar
%     axis equal off; drawnow
% 
% figure,contour(stream',60,'k-')
% 
% % title ('Streamline of Lid-Driven Cavity at Re 1000 with a cuboid within the centre')
% xlabel('x axis')
% ylabel('y axis')

% % Plot solution
%         nx = 10;
%         ny = 10;
%         xr = (0:1:Nx-1)*dx;
%         yr = (0:1:Ny-1)*dy;
%         [xc,yc] = ndgrid(xr,yr);
%         quiver(xc(1:nx:Nx, 1:ny:Ny), yc(1:nx:Nx, 1:ny:Ny), u(1:nx:Nx, 1:ny:Ny), v(1:nx:Nx, 1:ny:Ny),5);
%         axis([0 lenx 0 leny]);
%         drawnow  
