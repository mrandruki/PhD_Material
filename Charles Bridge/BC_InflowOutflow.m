 %    Specify the boundary conditions
 
 function [h,u,v]=BC_InflowOutflow(h,u,v,Q0,dy,dx,h0,u0,Lx,Ly)

% outflow
      v(Lx,:) = v(Lx-1, :); 
      u(Lx,:) = u(Lx-1,:);
      h(Lx,:) = h0; 
 
% inflow
%       for y=1:Ly
%           h(1,y)=h(2,y); % zero gradient
%       end
%       
%       for y = 1:Ly
%           v(1,y) = u0;
%       end
      
      for x=1:Lx
          h(x,1)=h(x,2); % zero gradient
      end
      
      for x = 1:Lx
          v(x,1) = u0;
          u(x,1) = u0;
      end

%       qi = 0; area = 0;
%       for y = 1:Ly-1
%          qi = qi+0.5*(h(1,y)+h(1,y+1))*dy*0.5*(u(1,y)+u(1,y+1));
%          area = area+0.5*(h(1,y)+h(1,y+1))*dy;
%       end
%       for y=1:Ly
%          u(1,y)= u(1,y)+(Q0-qi)/area;
%       end
      
      area = 0; qiv = 0; qiu = 0; 
      
      for x = 1:Lx-1
         qiv = qiv+0.5*(h(x,1)+h(x+1,1))*dx*0.5*(v(x,1)+v(x+1,1));
         qiu = qiu+0.5*(h(x,1)+h(x+1,1))*dx*0.5*(u(x,1)+u(x+1,1));
         area = area+0.5*(h(x,1)+h(x+1,1))*dx;
      end
      for x=1:Lx
         u(x,1) = u(x,1)+(Q0-qiu)/area; 
         v(x,1) = v(x,1)+(Q0-qiv)/area;
      end
     
% wall
      v(1,:) = 0.; 
      v(:,Ly) = 0;
      u(1,:) = 0;
  return