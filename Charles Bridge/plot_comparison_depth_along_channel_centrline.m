% Plot comparison depth along channel certrline

% The numerical data
load result.mat
plot (1:95,(h(1:95,51)-h0)/h0,'b')
hold on

% Theoretical data
%load ana.dat
ana=[1.54, 0.1861;		
1.67, 0.1894;	
1.71, 0.1916;	
1.725, 0.192;	
1.75, 0.1927;	
1.78, 0.1943;	
1.807, 0.196;	
1.835, 0.1982;	
1.85, 0.1993;	
1.86, 0.2015;	
1.87, 0.2037;	
1.88, 0.2059;	
1.89, 0.2087;	
2.11, 0.1694;	
2.13, 0.1696;	
2.14, 0.1696;	
2.165, 0.1694;	
2.193, 0.1694;	
2.275, 0.1691;	
2.33, 0.1694;	
2.44, 0.1696;	
2.66, 0.17]
xa=ana(1:13,1)*50;
ya=(ana(1:13,2)-h0)/h0;
plot (xa,ya,'ro');
xa=ana(16:22,1)*50;
ya=(ana(16:22,2)-h0)/h0;
plot (xa,ya,'ro');

load result.mat
plot (107:200,(h(107:200,51)-h0)/h0,'b')

% The cylinder
x=[95,95];
y=[-0.3,0.3]; 
plot(x,y)
x=[107,107];
y=[-0.3,0.3];
plot(x,y)
x=[95,107];
y=[0.3,0.3]; 
plot(x,y)



% Defined axises
xlabel('X(m)'), ylabel('(h-h0)/h0')
axis([70 130 -0.3 0.6])
set(gca,'XTickLabel',str2num(get(gca,'XTickLabel'))/50);


title ('Flow around cylinder: Comparison of the depth along channel centerline')
legend 'Lattice Boltzmann solution' 'Theoretical discharge'
hold off