% Determine solid values

function solid=trapezium_values(Lx,Ly,dx,dy,solid,trapezium)

pw = trapezium(1,:);
px = trapezium(2,:);
py = trapezium(3,:);
pz = trapezium(4,:);

for y=1:Ly
    for x=1:Lx
        if (y * dy) <= pw(2) && (y * dy) >= px(2)
            if (x * dx) >= pw(1) && (x * dx) >= px(1) && (x * dx) <= py(1) && (x * dx) <= pz(1)
                solid(x,y) = 2;
            else
                borderLeft = px(1) + (((y * dy) - px(2)) * ((pw(1) - px(1))/(pw(2) - px(2))));
                borderRight = pz(1) + (((y * dy) - pz(2)) * ((py(1) - pz(1))/(py(2) - pz(2))));

                if (x *dx) >= borderLeft && (x *dx) <= borderRight
                    solid(x,y) = 2;
                end
            end
        end
    end
end

return