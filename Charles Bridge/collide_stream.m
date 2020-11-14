% collide_stream(Xs-start,Xe-end,Ys-start,Ye-end)

function ftemp=collide_stream(Xs,Xe,Ys,Ye,tau,gacl,dt,Sx,Sy,u,v,e,ex,ey,feq,f,h,fb,solid)

for y=Ys:Ye
    yp=y+1;
    yn=y-1;
    for x=Xs:Xe
        if solid(x,y)<=1
           xp=x+1;
           xn=x-1;      
           if xp<=Xe &&solid(xp,y)<=1 
               ftemp(1,xp,y) = f(1,x,y)-(f(1,x,y)-feq(1,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(xp,y))*(ex(1)*Sx(x,y)+ey(1)*Sy(x,y))- dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(1)*u(x,y)+ey(1)*v(x,y));
           end
           if xp <= Xe && yp <= Ye &&solid(xp,yp)<=1
               ftemp(2,xp,yp) = f(2,x,y)-(f(2,x,y)-feq(2,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(xp,yp))*(ex(2)*Sx(x,y)+ey(2)*Sy(x,y))-dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(2)*u(x,y)+ey(2)*v(x,y));
           end
           if yp <= Ye &&solid(x,yp)<=1
               ftemp(3,x,yp) = f(3,x,y)-(f(3,x,y)-feq(3,x,y))/tau-dt/(6*e*e)*gacl*0.5*(h(x,y)+h(x,yp))*(ex(3)*Sx(x,y)+ey(3)*Sy(x,y))-dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(3)*u(x,y)+ey(3)*v(x,y));
           end
           if xn >= Xs && yp <= Ye && solid(xn,yp)<=1
               ftemp(4,xn,yp) = f(4,x,y)-(f(4,x,y)-feq(4,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(xn,yp))*(ex(4)*Sx(xn,y)+ey(4)*Sy(x,y))- dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(4)*u(x,y)+ey(4)*v(x,y));
           end
           if xn >= Xs && solid(xn,y)<=1 
               ftemp(5,xn,y) = f(5,x,y)-(f(5,x,y)-feq(5,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(xn,y))*(ex(5)*Sx(xn,y)+ey(5)*Sy(x,y))- dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(5)*u(x,y)+ey(5)*v(x,y));
           end
           if xn >= Xs && yn >= Ys && solid(xn,yn)<=1
               ftemp(6,xn,yn) = f(6,x,y)-(f(6,x,y)-feq(6,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(xn,yn))*(ex(6)*Sx(xn,y)+ey(6)*Sy(x,yn))- dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(6)*u(x,y)+ey(6)*v(x,y));
           end
           if yn >= Ys &&solid (x,yn)<=1
               ftemp(7,x,yn) = f(7,x,y) - (f(7,x,y)-feq(7,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(x,yn))*(ex(7)*Sx(x,y)+ey(7)*Sy(x,yn))- dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(7)*u(x,y)+ey(7)*v(x,y));
           end
           if xp <= Xe && yn >= Ys && solid(xp,yn)<=1
               ftemp(8,xp,yn) = f(8,x,y) - (f(8,x,y)-feq(8,x,y))/tau- dt/(6*e*e)*gacl*0.5*(h(x,y)+h(xp,yn))*(ex(8)*Sx(x,y)+ey(8)*Sy(x,yn))- dt/(6*e*e)*fb*sqrt(u(x,y)^2+v(x,y)^2)*(ex(8)*u(x,y)+ey(8)*v(x,y));
           end
           ftemp(9,x,y) = f(9,x,y) - (f(9,x,y)-feq(9,x,y))/tau;

        else
           ftemp(:,x,y) =0;
        end        
    end
end

return