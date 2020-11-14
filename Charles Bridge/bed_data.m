% This program is for bed data

function zb=bed_data(Lx,Ly,dx,bed_shape)

zb = zeros(Lx,Ly);
s = size(bed_shape);

for bs = 1:s-1
    y1 = bed_shape(bs,1);
    h1 = bed_shape(bs,2);
    
    y2 = bed_shape(bs+1,1);
    h2 = bed_shape(bs+1,2);
    
    if y2 > Ly
        y2 = Ly;
    end
    
    m = (h2-h1)/(y2-y1);
    
    for y = y1:y2
        for x = 1:Lx
            zb(x,y) = (m*(y-y1)) + h1;
        end
    end
end

if bed_shape(s,1)/dx < Ly
    for y = bed_shape(s,1)/dx:Ly
        for x = 1:Lx
            zb(x,y) = bed_shape(s,2);
        end
    end
end

return

% % sx=diff(zb,'x');
% % % 
% number_of_waves = 1000;
% factor = (360*number_of_waves*pi)/(Lx*180);  
% 
% for y=1:Ly
%     for x=1:Lx
%          zb(x,y) = 0.5 * sin((x-1)*dx*factor); 
% %                   sx*(x-1)*dx+0.0025; % negative
%       
%     end
% end         
            