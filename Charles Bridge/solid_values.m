% Determine solid values

function solid=solid_values(Lx,Ly,dx,dy,spb)

solid = zeros(Lx,Ly);

for sp = 1:size(spb)
    for y=1:Ly
        for x=1:Lx
            if spb(sp,1) == 1
                if sqrt(((x-1)*dx-spb(sp,2))^2+((y-1)*dy-spb(sp,3))^2 ) <= spb(sp,4)
                    solid(x,y)= 2;
                end
            elseif spb(sp,1) == 2
                if abs((x-1)*dx-spb(sp,2)) < spb(sp,4) && abs((y-1)*dy-spb(sp,3)) < spb(sp,4)
                    solid(x,y)= 2;
                end
            end
        end
    end
end

solid(1:Lx,1)=2;
solid(1:Lx,Ly)=2;

% solid(1:Lx,1)=1;
% solid(1:Lx,Ly)=1;

return