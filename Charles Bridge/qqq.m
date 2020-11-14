for i=1:Lx
QQ(i)=0
for j=2:Ly-1
QQ(i)=QQ(i)+0.5.*(h(1,j)+h(1,j+1)).*dy.*0.5.*(u(1,j)+u(1,j+1))
end
end

