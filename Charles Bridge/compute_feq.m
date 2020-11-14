%compute the equilibrium distribution function feq

function feq = compute_feq(Lx,Ly,h,u,v,e,ex,ey,gacl,solid)
for y=1:Ly
    for x=1:Lx
        if solid(x,y)<=1
            for a=1:8
                if mod(a,2)==0
                    feq(a,x,y)=gacl*h(x,y)^2/(24*e^2)+h(x,y)/(12*e^2)*(ex(a)*u(x,y)+ey(a)*v(x,y))+h(x,y)/(8*e^4)*(ex(a)*u(x,y)*ex(a)*u(x,y)+2*ex(a)*u(x,y)*ey(a)*v(x,y)+ey(a)*v(x,y)*ey(a)*v(x,y))-h(x,y)/(24*e^2)*(u(x,y)^2+v(x,y)^2);
                else
                    feq(a,x,y)=gacl*h(x,y)^2/(6*e^2)+h(x,y)/(3*e^2)*(ex(a)*u(x,y)+ey(a)*v(x,y))+h(x,y)/(2*e^4)*(ex(a)*u(x,y)*ex(a)*u(x,y)+2*ex(a)*u(x,y)*ey(a)*v(x,y)+ey(a)*v(x,y)*ey(a)*v(x,y))-h(x,y)/(6*e^2)*(u(x,y)^2+v(x,y)^2);
                end 
            end 
            feq(9,x,y)=h(x,y)-5*gacl*h(x,y)^2/(6*e^2)-2*h(x,y)/(3*e^2)*(u(x,y)^2+v(x,y)^2);  
        else
            feq(:,x,y)=0; 
        end
    end
end     

return
