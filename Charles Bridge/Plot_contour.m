% Plot contour

% The numerical data
load result.mat
Q=zeros(Lx,Ly);
for i=1:Lx
    for j=2:Ly-1
        Q(i,j)=(u(i,j)+u(i,j-1))*(h(i,j)+h(i,j-1))*dy/3;
    end
end
figure,contour(Q')
hold on

% The cuboids
x1 = [51 61 71 81 71 61 51];
y1 = [16.5 18.5 18.5 16.5 14.5 14.5 16.5];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [32.5 34.5 34.5 32.5 30.5 30.5 32.5];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [49 51 51 49 47 47 49];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [64 66 66 64 62 62 64];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [80.5 82.5 82.5 80.5 78.5 78.5 80.5];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [96.5 98.5 98.5 96.5 94.5 94.5 96.5];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [112.5 114.5 114.5 112.5 110.5 110.5 112.5];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [128 130 130 128 126 126 128];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

x1 = [51 61 71 81 71 61 51];
y1 = [144 146 146 144 142 142 144];
plot(x1,y1,'Color', 'r','LineWidth', 1.5);
hold on;

% Defined axises
xlabel('X(dm)'), ylabel('Y(dm)')
axis([80 201 0 80])
set(gca,'XTickLabel',str2num(get(gca,'XTickLabel'))/10);
set(gca,'YTickLabel',str2num(get(gca,'YTickLabel'))/10);

%title ('Flow around cylinder: contours of the water depth')
legend 'Contour'
hold off

surf(z)