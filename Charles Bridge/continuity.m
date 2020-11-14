% continuity test

function [Q]=continuity(Xs,Xe,Ys,Ye,h,u,dy)

for x=Xs:Xe
    Q(x)=0;
    for y=Ys:Ye-1
        Q(x)=Q(x)+0.5*(u(x,y)+u(x,y+1))*dy*0.5*(h(x,y)+h(x,y+1));
    end
end


return
