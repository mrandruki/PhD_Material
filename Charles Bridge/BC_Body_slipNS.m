% BC_Body (Xs,Xe,Lx,Ly)
% slip

function ftemp=BC_Body_slipNS(Xs,Xe,Ly, ftemp)

	  y = 1;
      for x=Xs:Xe
          ftemp(2,:,y) = ftemp(8,:,y); ftemp(3,:,y) = ftemp(7,:,y); ftemp(4,:,y) = ftemp(6,:,y);
      end

	  y = Ly;
      for x=Xs:Xe
          ftemp(8,:,y) = ftemp(2,:,y); ftemp(7,:,y) = ftemp(3,:,y); ftemp(6,:,y) = ftemp(4,:,y);
      end

return
