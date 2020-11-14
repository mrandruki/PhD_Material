! MRT-RLBM for N-S equations 

module ns2d

  implicit none

  integer:: max_t,x,y,a, Lx,Ly
  real:: dx,dy,ex0,ey0,nu_0,tau_0,dt,rho0,uo,fw, y0
  real:: sr = 1., c_smago = 0.1 !default values
  integer, allocatable, dimension(:,:):: solid
  real, allocatable, dimension(:):: ex,ey
  real, allocatable, dimension(:,:):: u,v,rho,Fx,Fy,tau, Mt, Mt_inv, S, MinvS
  real, allocatable, dimension(:,:,:):: f,feq,ftemp, omega, m, meq
  logical:: Tur_sgs = .false.

contains

  subroutine setup

    ! coordinate and other indices
    integer:: i
    real:: two_pi

    max_t = 1
    two_pi = 8.*atan(1.)

    ! compute allowed particle velocities
    ! neighbors are stored in the arrays as follows:

 
    ex0 = dx/dt
    ey0 = dy/dt

    !molecular viscosity
    nu_0 = ex0*ey0*dt*(2.*tau_0-1)/6.

    tau = tau_0

    ex(1) = ex0; ey(1) = 0.
    ex(2) = ex0; ey(2) = ey0
    ex(3) = 0.; ey(3) = ey0
    ex(4) = -ex0; ey(4) = ey0
    ex(5) = -ex0; ey(5) = 0.
    ex(6) = -ex0; ey(6) = -ey0
    ex(7) = 0.; ey(7) = -ey0
    ex(8) = ex0; ey(8) = -ey0
    ex(9) = 0.; ey(9) = 0.

    ! compute the equilibrium distribution function feq
    call compute_feq

    ! Set the initial distribution function to feq
    do y = 1, Ly
       do x = 1, Lx
          do a = 1, 9
             f(a,x,y) = feq(a,x,y)
          end do
       end do
    end do

  end subroutine setup


  subroutine collide_stream (Xs,Xe,Ys,Ye)

    ! coordinate and other indices
    integer:: a,b,xp,yp,xn,yn,Xs,Xe,Ys,Ye,ns,np,ic,jc
    real:: q_mft, tmp, ycor

    do y = Ys, Ye
       yp = y + 1
       yn = y - 1
      
       ycor = (y-1.)*dy - y0

       do x = Xs, Xe

          if (solid(x,y) == 2) cycle

!          m(:,x,y) = MATMUL(Mt, f(:,x,y))

!          m(:,x,y)=Mt*f(:,x,y); meq(:,x,y)=Mt*feq(:,x,y);

          m(:,x,y) = MATMUL(Mt, f(:,x,y))
          meq(:,x,y) = MATMUL(Mt, feq(:,x,y))

!         omega(:,x,y)=MinvS*(m(:,x,y)-meq(:,x,y));   
    
          omega(:,x,y) = MATMUL(MinvS, (m(:,x,y)-meq(:,x,y)))

          if (Tur_sgs) then

          ! tau is different,and must be computed on each point
          ! compute \sum_{i j=1)^b{\PI_{i j}^2}
          ! where \Pi_{i j}=sum_{\alpha=1}^D{e_{a,i} e_{a,j}(f_a-f^{eq}_a)}
          ! tmp is the momentum flux tensor \PI_{i j}

             q_mft = 0.
             do a = 1, 8
                tmp = ex(a)**2*(f(a,x,y)-feq(a,x,y))  ! \PI_xx
             end do
             q_mft = q_mft + tmp**2

             do a = 1, 8
                tmp = ex(a)*ey(a)*(f(a,x,y)-feq(a,x,y)) ! \PI_xy
             end do
             q_mft = q_mft + 2.*tmp**2

             do a = 1, 8
                tmp = ey(a)**2*(f(a,x,y)-feq(a,x,y)) ! \PI_yy
             end do
             q_mft = q_mft + tmp**2

           ! tmp=(sqrt(nu_0**2+18*c_smago**2*c_s2*sqrt(q_mft))-nu_0)/(6*c_smago**2*c_s2)

   tau(x,y) = (sqrt(tau_0**2+18*c_smago**2*sqrt(q_mft)/(rho(x,y)*ex0*ey0))+tau_0)/2.

          end if

          xp = x + 1
          xn = x - 1
          if (xp <= Xe) &
               & ftemp(1,xp,y) = f(1,x,y) - omega(1,x,y) &
               & + dt/(6.*ex0*ex0)*(ex(1)*Fx(x,y)+ey(1)*Fy(x,y))

          if (xp <= Xe .AND. yp <= Ye) &
               & ftemp(2,xp,yp) = f(2,x,y)  - omega(2,x,y) &
               & + dt/6.*(Fx(x,y)/ex(2)+Fy(x,y)/ey(2)) 

          if (yp <= Ye) &
               & ftemp(3,x,yp) = f(3,x,y) - omega(3,x,y) &
               & + dt/(6.*ey0*ey0)*(ex(3)*Fx(x,y)+ey(3)*Fy(x,y))

          if (xn >= Xs .AND. yp <= Ye) &
               & ftemp(4,xn,yp) = f(4,x,y) - omega(4,x,y) &
               & + dt/6.*(Fx(x,y)/ex(4)+Fy(x,y)/ey(4))

          if (xn >= Xs) &
               & ftemp(5,xn,y) = f(5,x,y)  - omega(5,x,y)&
               & + dt/(6.*ex0*ex0)*(ex(5)*Fx(x,y)+ey(5)*Fy(x,y))

          if (xn >= Xs .AND. yn >= Ys) &
               & ftemp(6,xn,yn) = f(6,x,y) - omega(6,x,y)&
               & + dt/6.*(Fx(x,y)/ex(6)+Fy(x,y)/ey(6))

          if (yn >= Ys)&
               & ftemp(7,x,yn) = f(7,x,y) - omega(7,x,y)&
               & + dt/(6.*ey0*ey0)*(ex(7)*Fx(x,y)+ey(7)*Fy(x,y))

          if (xp <= Xe .AND. yn >= Ys)&
               & ftemp(8,xp,yn) = f(8,x,y)  - omega(8,x,y) &
               & + dt/6.*(Fx(x,y)/ex(8)+Fy(x,y)/ey(8)) 

          ftemp(9,x,y) = f(9,x,y) - omega(9,x,y) 

       end do
    end do  


    
    return
  end subroutine collide_stream


  subroutine BC_InflowOutflow (InOut,Xc,Ys,Ye)
    integer:: Xc,Ys,Ye 
    character:: InOut*3

    ! Specify the inflow/outflow BCs y=Ye (Ys < y < Ye)

    x = Xc

    if (InOut == "in") then
       do y = Ys, Ye   
          if (solid(x,y) == 2) cycle
          rho(x,y) = ex0/(ex0-u(x,y))*(ftemp(3,x,y)+ftemp(7,x,y)+ftemp(9,x,y)&
               & +2.*(ftemp(4,x,y)+ftemp(5,x,y)+ftemp(6,x,y)))
          ftemp(1,x,y) = ftemp(5,x,y)+2.*rho(x,y)*u(x,y)/(3.*ex0)
          ftemp(2,x,y) = rho(x,y)*u(x,y)/(6.*ex0)+ftemp(6,x,y)&
               & +(ftemp(7,x,y)-ftemp(3,x,y))/2.
          ftemp(8,x,y) = rho(x,y)*u(x,y)/(6.*ex0)+ftemp(4,x,y)&
               & +(ftemp(3,x,y)-ftemp(7,x,y))/2.
       end do
    else if (InOut == "out") then
       do y = Ys, Ye  
          if (solid(x,y) == 2) cycle 
          u(x,y) = -ex0+ex0/rho(x,y)*(ftemp(3,x,y)+ftemp(7,x,y)+ftemp(9,x,y)&
               & +2.*(ftemp(1,x,y)+ftemp(2,x,y)+ftemp(8,x,y)))
          ftemp(5,x,y) = ftemp(1,x,y)-2.*rho(x,y)*u(x,y)/(3.*ex0)
          ftemp(4,x,y) = -rho(x,y)*u(x,y)/(6.*ex0)+ftemp(8,x,y)&
               & +(ftemp(7,x,y)-ftemp(3,x,y))/2.
          ftemp(6,x,y) = -rho(x,y)*u(x,y)/(6.*ex0)+ftemp(2,x,y)&
               & +(ftemp(3,x,y)-ftemp(7,x,y))/2.
       end do
    end if
 
    return
  end subroutine BC_InflowOutflow
  


  subroutine BC_NorthSouth (Xs,Xe,Ys,Ye)

    integer:: Xs,Xe,Ys,Ye

    ! Specify the north BC y=Ye,Ys (Xs < x < Xe)

    y = Ye
    do x = Xs, Xe
       rho(x,y) = ftemp(1,x,y)+ftemp(5,x,y)+ftemp(9,x,y)&
            & +2.*(ftemp(2,x,y)+ftemp(3,x,y)+ftemp(4,x,y))
       ftemp(7,x,y) = ftemp(3,x,y)
       ftemp(6,x,y) = -rho(x,y)*u(x,y)/(2.*ex0)+ftemp(2,x,y)&
            & +0.5*(ftemp(1,x,y)-ftemp(5,x,y))
       ftemp(8,x,y) = rho(x,y)*u(x,y)/(2.*ex0)+ftemp(4,x,y)&
            & +0.5*(ftemp(5,x,y)-ftemp(1,x,y))
    end do
    y = Ys
    do x = Xs, Xe
       rho(x,y) = ftemp(1,x,y)+ftemp(5,x,y)+ftemp(9,x,y)&
            & +2.*(ftemp(6,x,y)+ftemp(7,x,y)+ftemp(8,x,y))
       ftemp(3,x,y) = ftemp(7,x,y)
       ftemp(2,x,y) = rho(x,y)*u(x,y)/(2.*ex0)+ftemp(6,x,y)&
            & +0.5*(ftemp(5,x,y)-ftemp(1,x,y))
       ftemp(4,x,y) = -rho(x,y)*u(x,y)/(2.*ex0)+ftemp(8,x,y)&
            & +0.5*(ftemp(1,x,y)-ftemp(5,x,y))
    end do

    return
  end subroutine BC_NorthSouth


  subroutine BC_West (Xw,Ys,Ye)

    integer:: Xw,Ys,Ye

    ! Specify the north BC y=Ye,Ys (Xs < x < Xe)

    x = Xw
    do y = Ys, Ye
       ftemp(1,x,y) = ftemp(5,x,y)
       ftemp(2,x,y) = ftemp(6,x,y)
       ftemp(8,x,y) = ftemp(4,x,y)
    end do

    return
  end subroutine BC_West

  subroutine BC_East (Xe,Ys,Ye)

    integer:: Xe,Ys,Ye

    ! Specify the north BC y=Ye,Ys (Xs < x < Xe)

    x = Xe
    do y = Ys, Ye
       ftemp(5,x,y) = ftemp(1,x,y)
       ftemp(6,x,y) = ftemp(2,x,y)
       ftemp(4,x,y) = ftemp(8,x,y)
    end do

    return
  end subroutine BC_East
  

  subroutine BC_CornerSW (xsw,ysw)

    integer:: xsw,ysw

    !Corner south west point
    x = xsw; y = ysw
    ftemp(1,x,y) = 2.*rho(x,y)*u(x,y)/(3.*ex0)+ftemp(5,x,y)
    ftemp(3,x,y) = ftemp(7,x,y)
    ftemp(2,x,y) = rho(x,y)*u(x,y)/(6.*ex0)+ftemp(6,x,y)
    ftemp(4,x,y) = 0.5*(rho(x,y)-(ftemp(1,x,y)+ftemp(2,x,y)+ftemp(3,x,y)&
         & +ftemp(5,x,y)+ftemp(6,x,y)+ftemp(7,x,y)+ftemp(9,x,y)))
    ftemp(8,x,y) = ftemp(4,x,y)

    return
  end subroutine BC_CornerSW
  

  subroutine BC_CornerNW (xnw,ynw)
    integer:: xnw,ynw

    !Corner northth west point
    x = xnw; y = ynw
    ftemp(1,x,y) = 2.*rho(x,y)*u(x,y)/(3.*ex0)+ftemp(5,x,y)
    ftemp(7,x,y) = ftemp(3,x,y)
    ftemp(8,x,y) = rho(x,y)*u(x,y)/(6.*ex0)+ftemp(4,x,y)
    ftemp(2,x,y) = 0.5*(rho(x,y)-(ftemp(1,x,y)+ftemp(3,x,y)+ftemp(4,x,y)&
         & +ftemp(5,x,y)+ftemp(7,x,y)+ftemp(8,x,y)+ftemp(9,x,y)))
    ftemp(6,x,y) = ftemp(2,x,y)

    return
  end subroutine BC_CornerNW
  
  subroutine BC_CornerSE (xse,yse)
    integer:: xse,yse

    !Corner South-Eastt point

    x = xse; y = yse
    ftemp(5,x,y) = -2.*rho(x,y)*u(x,y)/(3.*ex0)+ftemp(1,x,y)
    ftemp(3,x,y) = ftemp(7,x,y) 
    ftemp(4,x,y) = -rho(x,y)*u(x,y)/(6.*ex0)+ftemp(8,x,y)
    ftemp(2,x,y) = 0.5*(rho(x,y)-(ftemp(1,x,y)+ftemp(3,x,y)+ftemp(4,x,y)&
         & +ftemp(5,x,y)+ftemp(7,x,y)+ftemp(8,x,y)+ftemp(9,x,y)))
    ftemp(6,x,y) = ftemp(2,x,y)

    return
  end subroutine BC_CornerSE
  
  subroutine BC_CornerNE (xne,yne)
    integer:: xne,yne

    !Corner South-Eastt point
    x = xne; y = yne
    ftemp(5,x,y) = -2.*rho(x,y)*u(x,y)/(3.*ex0)+ftemp(1,x,y)
    ftemp(7,x,y) = ftemp(3,x,y)
    ftemp(6,x,y) = -rho(x,y)*u(x,y)/(6.*ex0)+ftemp(2,x,y)
    ftemp(4,x,y) = 0.5*(rho(x,y)-(ftemp(1,x,y)+ftemp(2,x,y)+ftemp(3,x,y)&
         & +ftemp(5,x,y)+ftemp(6,x,y)+ftemp(7,x,y)+ftemp(9,x,y)))
    ftemp(8,x,y) = ftemp(4,x,y)

  end subroutine BC_CornerNE


  subroutine No_Slipx (bound,Yc,Xs,Xe)

    integer:: Yc,Xs,Xe
    character:: bound*5

    y = Yc
    if (bound == "south") then 
       do x = Xs, Xe
          ftemp(2,x,y) = ftemp(6,x,y)
          ftemp(3,x,y) = ftemp(7,x,y)
          ftemp(4,x,y) = ftemp(8,x,y)
       end do
    end if

    if (bound == "north") then
       do x = Xs, Xe
          ftemp(6,x,y) = ftemp(2,x,y)
          ftemp(7,x,y) = ftemp(3,x,y)
          ftemp(8,x,y) = ftemp(4,x,y)
       end do
    end if


    return
  end subroutine No_Slipx

  
  subroutine Slipx (bound,Yc,Xs,Xe)

    integer:: Yc,Xs,Xe
    character:: bound*5

    y = Yc
    if (bound == "south") then 
       do x = Xs-1, Xe+1      
          ftemp(2,x,y) = ftemp(8,x,y)
          ftemp(3,x,y) = ftemp(7,x,y)
          ftemp(4,x,y) = ftemp(6,x,y)
       end do
    end if

    if (bound == "north") then
       do x = Xs-1, Xe+1          
          ftemp(6,x,y) = ftemp(4,x,y)
          ftemp(7,x,y) = ftemp(3,x,y)
          ftemp(8,x,y) = ftemp(2,x,y)
       end do
    end if
    
    return
  end subroutine Slipx



  subroutine Semi_Slipx (bound,Yc,Xs,Xe)

    integer:: Xs,Xe,Yc
    character:: bound*5

    y = Yc
    if (bound == "south") then 
       do x = Xs-1, Xe+1   
          do a = 6, 8
             ftemp(a,x,y) = ftemp(a,x,y)- dt/(6.*ex0*ex0)*fw* &
                  & sqrt(u(x,y)**2+v(x,y)**2)* &
                  & (ex(a)*u(x,y)+ey(a)*v(x,y)) 
          end do
          ftemp(2,x,y) = ftemp(8,x,y)
          ftemp(3,x,y) = ftemp(7,x,y)
          ftemp(4,x,y) = ftemp(6,x,y)
       end do
    end if

    if (bound == "north") then
       do x = Xs-1, Xe+1   
          do a = 2, 4
             ftemp(a,x,y) = ftemp(a,x,y)- dt/(6.*ex0*ex0)*fw* &
                  & sqrt(u(x,y)**2+v(x,y)**2)* &
                  & (ex(a)*u(x,y)+ey(a)*v(x,y)) 
          end do
          ftemp(6,x,y) = ftemp(4,x,y)
          ftemp(7,x,y) = ftemp(3,x,y)
          ftemp(8,x,y) = ftemp(2,x,y)
       end do

    end if

    return
  end subroutine Semi_Slipx



  subroutine BC_periodX (bound,Xc,Ys,Ye)

    integer:: Xc,Ys,Ye
    character:: bound*2


    x = Xc
    if (bound == "in") then
       do y = Ys, Ye
          ftemp(1,x,y) = ftemp(1,Lx,y)
          ftemp(2,x,y) = ftemp(2,Lx,y)
          ftemp(8,x,y) = ftemp(8,Lx,y)
       end do
    end if

    if (bound == "ou") then
       do y = Ys, Ye
          ftemp(4,x,y) = ftemp(4,1,y)
          ftemp(5,x,y) = ftemp(5,1,y)
          ftemp(6,x,y) = ftemp(6,1,y)
       end do
    end if

    return
  end subroutine BC_periodX


  


  subroutine solution (Xs,Xe,Ys,Ye)

    integer:: Xs,Xe,Ys,Ye

    !Set the global f
    do y = Ys, Ye
       do x = Xs, Xe
          if (solid(x,y) == 2) cycle
          do a = 1, 9
             f(a,x,y) = ftemp(a,x,y)
!!             if (x == Lx) f(a,x,y) = feq(a,x,y)
          end do
       end do
    end do

    !compute the velocity and depth
    call compute_h_and_u (Xs,Xe,Ys,Ye)

    return
  end subroutine solution


  subroutine compute_h_and_u (Xs,Xe,Ys,Ye)

    integer:: Xs,Xe,Ys,Ye

    do y = Ys, Ye
       do x = Xs, Xe
          if (solid(x,y) == 2) cycle
          rho(x,y) = 0.0
          u(x,y) = 0.0
          v(x,y) = 0.0
          do a = 1, 9
             rho(x,y) = rho(x,y) + f(a,x,y)
             u(x,y) = u(x,y) + ex(a)*f(a,x,y)
             v(x,y) = v(x,y) + ey(a)*f(a,x,y)
          end do
          u(x,y) = u(x,y)/rho(x,y)
          v(x,y) = v(x,y)/rho(x,y)
       end do
    end do

    return
  end subroutine compute_h_and_u


  subroutine compute_feq

  real:: beta

  beta = 1./6.

  beta = 1./9.

  beta = 1./12.

  beta = 1./15.  ! ok for dx=2dy

  beta = 1./18.

  beta = 1./30.

!  beta = 1./36.

!  beta = 1./72.

    do y = 1, Ly
       do x = 1, Lx
          if (solid(x,y) == 2) cycle
          do a = 1, 8

             if (mod(a,2) == 0) then
                feq(a,x,y) = rho(x,y)*( 1./12.*(u(x,y)/ex(a)+ v(x,y)/ey(a)) &
                     & +1./4.*u(x,y)*v(x,y)/(ex(a)*ey(a)) )
             else if (a == 1 .or. a == 5) then
                feq(a,x,y) = rho(x,y)*( beta*ey0/ex0 + &
                     & u(x,y)/(3*ex(a))+ u(x,y)*u(x,y)/(2*ex0*ex0) )
             else if (a == 3 .or. a == 7) then
                feq(a,x,y) = rho(x,y)*( beta*ex0/ey0 + &
                     & v(x,y)/(3*ey(a))+ v(x,y)*v(x,y)/(2*ey0*ey0) )
             end if

          end do

          feq(9,x,y) = rho(x,y)- rho(x,y)*( 2*beta*ey0/ex0 + u(x,y)*u(x,y)/(ex0*ex0) ) &
               &               - rho(x,y)*( 2*beta*ex0/ey0 + v(x,y)*v(x,y)/(ey0*ey0) )

       end do
    end do

    return
  end subroutine compute_feq
 
end module ns2d
