% Plot velocity vectors

% The numerical data
load result.mat
nskpx = 3;
nskpy = 1;
for x=1:nskpx:Lx
    xn=(x-1)*dx;  
    for y=1:nskpy:Ly
            yn=(y-1)*dy;
            quiver(xn,yn,u(x,y),v(x,y),'k','filled','LineWidth',1)
            hold on
    end
end

% The cylinder
r = 0.3;
sita=0:pi/20:2*pi;
plot(10+r'*cos(sita),4.75+r'*sin(sita),'LineWidth', 2);
hold on

r = 0.3;
sita=0:pi/20:2*pi;
plot(10+r'*cos(sita),2.75+r'*sin(sita),'LineWidth', 2);
hold on

r = 0.3;
sita=0:pi/20:2*pi;
plot(10+r'*cos(sita),0.75+r'*sin(sita),'LineWidth', 2);
hold on

 

% Set axises
xlabel('x(m)'), ylabel('y(m)')
axis([8 18 0 5.5])

title ('Flow around cylinder: velocity vectors')
legend 'Lattice Boltzmann solution' 
hold off