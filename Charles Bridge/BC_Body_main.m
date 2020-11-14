% BC_Body (Xs,Xe,Ys,Ye)
% no-slip

function ftemp=BC_Body_main(Xs,Xe,Ys,Ye,solid,ftemp)

   for y = Ys+1:Ye-1
       for x = Xs+1: Xe-1
          if (solid(x,y) == 1) || (solid(x,y) == 0)
              if (solid(x+1,y) == 2) 
                  ftemp(5,x,y) = ftemp(1,x,y);
              end
              if (solid(x+1,y+1) == 2)
                  ftemp(6,x,y) = ftemp(2,x,y);
              end
              if (solid(x,y+1) == 2) 
                  ftemp(7,x,y) = ftemp(3,x,y);
              end
              if (solid(x-1,y+1) == 2) 
                  ftemp(8,x,y) = ftemp(4,x,y);
              end
              if (solid(x-1,y) == 2)
                  ftemp(1,x,y) = ftemp(5,x,y);
              end
              if (solid(x-1,y-1) == 2)
                  ftemp(2,x,y) = ftemp(6,x,y);
              end
              if (solid(x,y-1) == 2)
                  ftemp(3,x,y) = ftemp(7,x,y);
              end
              if (solid(x+1,y-1) == 2)
                  ftemp(4,x,y) = ftemp(8,x,y);          
              end
          end
       end
    end

return
