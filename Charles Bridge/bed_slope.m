% Bed Slop information

function [Sx,Sy]=bed_slope(zb,Lx,Ly,dx,dy)

% S-Slope
for x=1:Lx-1;
    for y=1:Ly-1;
        Sx(x,y)=0;
        Sy(x,y)=0;
    end
end

% Deal with boundary slopes
for i=1:Ly-1
    Sx(Lx,i)=0;
    Sy(Lx,i)=0;
end

for j=1:Lx-1
    Sx(j,Ly)=0;
    Sy(j,Ly)=0;
end

Sx(Lx,Ly)=0;
Sy(Lx,Ly)=0;

return