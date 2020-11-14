 % Specify the inlet and outlet boundary conditions
 
 function f=BC_InOut(Lx,Ly,f)
 
	f(:,1,:) = f(:,2,:);      % important
	f(:,Lx,:) = f(:,Lx-1,:);  % important
    f(:,:,1) = f(:,:,2);
        
  return