% Solution (Xs,Xe,Ys,Ye)

function [h,u,v,f]=solution(Xs,Xe,Ys,Ye,solid,ftemp,ex,ey,f)

for y = Ys:Ye
    for x = Xs: Xe
          for a = 1: 9
             f(a,x,y) = ftemp(a,x,y);
          end
    end
end

% Compute the velocity and depth
for y = Ys:Ye
    for x = Xs:Xe
        if solid(x,y)<=1
          h(x,y) = 0;
          sum_u(x,y) = 0;
          sum_v(x,y) = 0;
          for a = 1:9
             h(x,y) = h(x,y) + f(a,x,y);
             sum_u(x,y) = sum_u(x,y) + ex(a)*f(a,x,y);
             sum_v(x,y) = sum_v(x,y) + ey(a)*f(a,x,y);
          end
          u(x,y) = sum_u(x,y)/h(x,y);
          v(x,y) = sum_v(x,y)/h(x,y); 
        elseif solid(x,y)==2
          h(x,y) = 0;v(x,y)=0;  u(x,y)=0;  
        end
    end
end

return