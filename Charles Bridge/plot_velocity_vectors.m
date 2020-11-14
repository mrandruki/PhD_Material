% Plot velocity vectors

% The numerical data
load result.mat

n = 3;
xr=(0:1:Lx-1)*dx;
yr=(0:1:Ly-1)*dy;
[xc,yc] = ndgrid(xr,yr);
quiver(xc(1:n:Lx, 1:1:Ly),yc(1:n:Lx, 1:1:Ly),u(1:n:Lx, 1:1:Ly),v(1:n:Lx, 1:1:Ly), 3);
hold on

% The cuboids
x1 = [5 6 7 8 7 6 5];
y1 = [1.6 1.8 1.8 1.6 1.4 1.4 1.6];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [3.2 3.4 3.4 3.2 3 3 3.2];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [4.8 5 5 4.8 4.6 4.6 4.8];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [6.4 6.6 6.6 6.4 6.2 6.2 6.4];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [8 8.2 8.2 8 7.8 7.8 8];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [9.6 9.8 9.8 9.6 9.4 9.4 9.6];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [11.2 11.4 11.4 11.2 11 11  11.2];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [12.8 13 13 12.8 12.6 12.6 12.8];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [5 6 7 8 7 6 5];
y1 = [14.4 14.6 14.6 14.4 14.2 14.2 14.4];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

% Set axises
xlabel('x(m)'), ylabel('y(m)')
axis([-1 21 -1 9])
pbaspect([2.2 1 1])

title ('Charles Bridge, Prague, Czech Republic')
legend 'Lattice Boltzmann solution' 
hold off