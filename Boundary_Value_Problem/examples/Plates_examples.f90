module Plates_examples 
    
    use Boundary_value_problems 
    use Finite_differences
    
    implicit none
    
    contains  
    
    
    
!*******************************************************************************************************************************************
! system Linear_2D
!*******************************************************************************************************************************************
subroutine Test_Non_Linear_system_BVP2D

    integer, parameter :: Nx = 20, Ny = 20, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    integer :: i, j
    real :: a=0, b=1
    real :: PI = 4 * atan(1.0)

    x = [ (a + (b-a)*i/Nx, i=0, Nx) ]
    y = [ (a + (b-a)*j/Ny, j=0, Ny) ]
    U = 1      

    !call Non_Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 11, N_variables = 2, & 
    !                                        Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          
   
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \nabla^4 u =   50  \quad sin(2 \pi x) \quad sin( 2 \pi y )  $ ", 2)
    CALL TITLIN (" $ u(0,y) = 0, \quad u(1,y) = 0, \quad u_y(x,0) = 0, \quad u(x,1) = 0 $  ", 4)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)

contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
        
       real :: v, vxx, vyy, w, wxx, wyy
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) 
        w = u(2); wxx = uxx(2); wyy = uyy(2)
        
       
        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - 50 * sin(2*PI*x) * sin(2*PI*y)
        
end function

!********* Boundary conditions *********
function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
        elseif (x==b) then
             BCs(1) = u(1)
             BCs(2) = u(2)
        elseif (y==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
        elseif (y==b) then
            BCs(1) = u(1)
            BCs(2) = u(2)
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " a, b=", a, b
            stop 
        endif

    end function

end subroutine 




!*******************************************************************************************************************************************
! system Linear_2D
!*******************************************************************************************************************************************
subroutine Biharmonic2D

    integer, parameter :: Nx = 20, Ny = 20, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    integer :: i, j
    real :: a=0, b=1, q0 = 1
    real :: PI = 4 * atan(1.0)

    x = [ (a + (b-a)*i/Nx, i=0, Nx) ]
    y = [ (a + (b-a)*j/Ny, j=0, Ny) ]
    U = 1      

    call Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 11, N_variables = 2, & 
                                 Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          

contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
        
       real :: v, vxx, vyy, w, wxx, wyy
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) 
        w = u(2); wxx = uxx(2); wyy = uyy(2)
        
       
        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q0
        
end function

!********* Boundary conditions *********
function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
        elseif (x==b) then
             BCs(1) = u(1)
             BCs(2) = u(2)
        elseif (y==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
        elseif (y==b) then
            BCs(1) = u(1)
            BCs(2) = u(2)
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " a, b=", a, b
            stop 
        endif

    end function

end subroutine 

!*******************************************************************************************************************************************
! system Linear_3D
!*******************************************************************************************************************************************
subroutine Test_Linear_system_BVP3D

    integer, parameter :: Nx = 10, Ny = 10, Nz = 10,  Nv = 2 
    real :: x(0:Nx), y(0:Ny), z(0:Nz), U(0:Nx, 0:Ny, 0:Nz, Nv)
    integer :: i, j, k 
    real :: a=0, b=1
    real :: PI = 4 * atan(1.0)

    x = [ (a + (b-a)*i/Nx, i=0, Nx) ]
    y = [ (a + (b-a)*j/Ny, j=0, Ny) ]
    z = [ (a + (b-a)*k/Nz, k=0, Nz) ]
    U = 1      

    !call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, z_nodes = z, Order = 4, N_variables = 2, & 
    !                                    Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          
   
    write(*,*) " maxval U =", maxval( U(:,:,:,1) )     
    write(*,*) " minval U =", minval( U(:,:,:,1) ) 
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \nabla^4 u =   50  \quad sin(2 \pi x) \quad sin( 2 \pi y ) \quad sin( 2 \pi z )  $ ", 2)
    CALL TITLIN (" $ u(0,y,z) = 0, \quad u(1,y,z) = 0, \quad u_y(x,0,z) = 0, \quad u(x,1,z) = 0 $  ", 4)
    
    call qplcon( U(:,:, 0, 1), Nx+1, Ny+1, 20)
    call qplcon( U(:,:, 2, 1), Nx+1, Ny+1, 20)
    call qplcon( U(:,:, Nz, 1), Nx+1, Ny+1, 20)
    

contains
  
!********* Function *********
function L(x, y, z, u, ux, uy, uz, uxx, uyy, uzz, uxy, uxz, uyz)
  real, intent(in) :: x, y, z, u(:), ux(:), uy(:), uz(:), uxx(:), uyy(:), uzz(:), uxy(:), uxz(:), uyz(:)
  real :: L(size(u)) 
        
       real :: v, vxx, vyy, vzz, w, wxx, wyy, wzz
       
        v = u(1); vxx = uxx(1); vyy = uyy(1); vzz = uzz(1)     
        w = u(2); wxx = uxx(2); wyy = uyy(2); wzz = uzz(2) 
        
       
        L(1) =  vxx + vyy  + vzz  -50 * sin(2*PI*x) * sin(2*PI*y)  * sin(2*PI*z) ! ! 
        L(2) =  wxx + wyy + wzz - 50 * sin(2*PI*x) * sin(2*PI*y) * sin(2*PI*z) 
        
end function   

!********* Boundary conditions *********
function BCs(x, y, z, u, ux, uy, uz)
    
        real, intent(in) :: x, y, z, u(:), ux(:), uy(:), uz(:) 
        real :: BCs(size(u))

        if (x==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==b) then
             BCs(1) = u(1)
             BCs(2) = u(2)
             
        elseif (y==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (y==b) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (z==a) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (z==b) then
            BCs(1) = u(1)
            BCs(2) = u(2)    
            
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " a, b=", a, b
            stop 
        endif

    end function

end subroutine 



!*******************************************************************************************************************************************
! Non_Linear_2D
!*******************************************************************************************************************************************
   

subroutine Test_Non_Linear_BVP2D

    integer, parameter :: Nx = 20, Ny = 20
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny)
    integer :: i, j
    real :: a=0, b=1
    real :: pi = 4 * atan(1.0)

    x = [ (a + (b-a)*i/Nx, i=0, Nx) ]
    y = [ (a + (b-a)*j/Ny, j=0, Ny) ]
    U = 1      

    call Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 5,  & 
                                 Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          
   
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ ( u_{xx} + u_{yy} ) u = 0  $ ", 2)
    CALL TITLIN (" $ u(0,y) = 0, \quad u(1,y) = y, \quad u_y(x,0) = 0, \quad u(x,1) = x $  ", 4)
    call QPLClr( U, Nx+1, Ny+1)
       
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ ( u_{xx} + u_{yy} ) u = 0  $ ", 2)
    CALL TITLIN (" $ u(0,y) = 0, \quad u(1,y) = y, \quad u_y(x,0) = 0, \quad u(x,1) = x $  ", 4)
    call qplcon( U, Nx+1, Ny+1, 20)

contains

!********* Function *********

    real function L(x, y, u, ux, uy, uxx, uyy, uxy)
    
        real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy

        L = ( uxx + uyy) * u

    end function

!********* Bounds *********

    real function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u, ux, uy

        if (x==a) then
            BCs = u
        elseif (x==b) then
            BCs = u - y
        elseif (y==a) then
            BCs = uy
        elseif (y==b) then
            BCs = u - x
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " a, b=", a, b
            stop 
        endif

    end function

end subroutine 









!*******************************************************************************************************************************************
! system Linear_3D
!*******************************************************************************************************************************************
subroutine Test_plate3D

    integer, parameter :: Nx = 15, Ny = 15, Nz = 2,  Nv = 3 
    real :: x(0:Nx), y(0:Ny), z(0:Nz), U(0:Nx, 0:Ny, 0:Nz, Nv)
    integer :: i, j, k 
    real :: x0, xf, y0, yf, z0, zf
    real :: Lx = 1.3, Ly = 3, Lz = 1.6d-2
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_pillar = 0 , y_pillar = 0.50 , R_pillar = 20d-3 

!************* Uniform grid ********************************   
    x0 = -0.5*Lx ;  y0 = -0.5*Ly ; z0 = -0.5*Lz
    xf =  0.5*Lx ;  yf =  0.5*Ly ; zf =  0.5*Lz
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]
    z = [ (z0 + (zf-z0)*k/Nz, k=0, Nz) ]

!***********************************************************    
    mu = 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    K = (mu + lambda)/mu  
    U = 1      

    !call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, z_nodes = z, Order = 4, N_variables = Nv, & 
    !                                    Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          
   
    write(*,*) " maxval U =", maxval( U(:,:,:,3) )     
    write(*,*) " minval U =", minval( U(:,:,:,3) ) 
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $    $ ", 2)
    CALL TITLIN (" $ u(0,y,z) = 0, \quad u(1,y,z) = 0, \quad u_y(x,0,z) = 0, \quad u(x,1,z) = 0 $  ", 4)
    
    call qplcon( U(:,:,  0, 3), Nx+1, Ny+1, 20)
    call qplcon( U(:,:,  Nz/2, 3), Nx+1, Ny+1, 20)
    call qplcon( U(:,:,  Nz, 3), Nx+1, Ny+1, 20)
    
    call qplcon( U(:, Ny/2,:,3), Nx+1, Nz +1,20)
    

contains
  
!********* Function *********
function L(x, y, z, u, ux, uy, uz, uxx, uyy, uzz, uxy, uxz, uyz)
  real, intent(in) :: x, y, z, u(:), ux(:), uy(:), uz(:), uxx(:), uyy(:), uzz(:), uxy(:), uxz(:), uyz(:)
  real :: L(size(u)) 
  real :: P, Pxx, Pyy, Pzz, Pxy, Pxz, Pyz
  real :: Q, Qxx, Qyy, Qzz, Qxy, Qxz, Qyz
  real :: R, Rxx, Ryy, Rzz, Rxy, Rxz, Ryz



    P = u(1)  ; Pxx = uxx(1); Pyy = uyy(1); Pzz = uzz(1)  
  Pxy = uxy(1); Pxz = uxz(1); Pyz = uyz(1)   
 
    Q = u(2)  ; Qxx = uxx(2); Qyy = uyy(2); Qzz = uzz(2)     
  Qxy = uxy(2); Qxz = uxz(2); Qyz = uyz(2)   

    R = u(3)  ; Rxx = uxx(3); Ryy = uyy(3); Rzz = uzz(3)     
  Rxy = uxy(3); Rxz = uxz(3); Ryz = uyz(3)   
 
        
  L(1) = ( Pxx + Pyy + Pzz ) + K*( Pxx + Qxy + Rxz ) 
  L(2) = ( Qxx + Qyy + Qzz ) + K*( Pxy + Qyy + Ryz ) 
  L(3) = ( Rxx + Ryy + Rzz ) + K*( Pxz + Qyz + Rzz ) 
        
end function   

!********* Boundary conditions *********
function BCs(x, y, z, u, ux, uy, uz)
    
        real, intent(in) :: x, y, z, u(:), ux(:), uy(:), uz(:) 
        real :: BCs(size(u))
!~         real :: pillar_region
        
!~   pillar_region = (x - x_pillar)**2 +(y - y_pillar)**2 - R_pillar**2

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            
        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            
        elseif (z==z0) then
            BCs(1) = uz(1) + ux(3)
            BCs(2) = uz(2) + uy(3)
            BCs(3) = (lambda + 2*mu)* ux(1) + lambda*( uy(2) + uz(3) ) + p_atm
           
        elseif (z==zf) then
        
!~           if (pillar_region > 0) then
!~             BCs(1) = uz(1) + ux(3)
!~             BCs(2) = uz(2) + uy(3)
!~             BCs(3) = (lambda + 2*mu)* ux(1) + lambda*( uy(2) + uz(3) ) + p_atm - rho_water * g *y 
            
!~           else
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            
!~           endif 
                
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            stop 
        endif

    end function

end subroutine 








!*******************************************************************************************************************************************
! Pillars system Navier 2D
!*******************************************************************************************************************************************


subroutine Test_Pillars_Navier_system

    integer, parameter :: Nx = 36, Ny = 36, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 8d-3, Lx = 1.3, Ly = 3
    real :: x0 , xf , y0 , yf
    integer :: i, j, i_pillar, j_pillar
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_pillar, y_pillar, R_pillar
    real :: dmin 
    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly
    x_pillar = Lx/2 
    y_pillar = 2*Ly/3 
    dmin = max( Lx/Nx, Ly/Ny) 
    R_pillar = max( 0.2, dmin ) 
    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

  !  call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 4, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
     
  call Grid_Initialization( "nonuniform", "x", 4, x )
  call Grid_Initialization( "nonuniform", "y", 4, y )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      
  open(62, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Pillars_Navier_Cartesian\Pillars_Navier_Cartesian.plt')
write (62,*) 'VARIABLES = "X", "Y", "Desplazamiento (m)" , "Sxx", "Syy", "Sxy"'
write (62,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
do  j= 0, Ny
  do i = 0, Nx
 write (62,*) x(i), y(j), U(i,j,1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
end do
close (62)

 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Pillars_Navier_Cartesian\Pillars_Navier_Cartesian2.plt')
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy
  real :: pillar_region , q
    
  pillar_region = (x - x_pillar)**2 +(y - y_pillar)**2 - R_pillar**2
  q = ( - rho_water * g * (y - y0))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2) 
       
    if (pillar_region > 0) then   
        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
       
     else
        L(1) = vxx +vyy - w 
        L(2) = v 
      
    endif  
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 




!*******************************************************************************************************************************************
! Pillars Non Linear system Navier 2D
!*******************************************************************************************************************************************


subroutine Test_Pillars_Non_Linear_Navier_system

    integer, parameter :: Nx = 36, Ny = 36, Nv = 4
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    real :: Ux(0:Nx, 0:Ny, Nv), Uy(0:Nx, 0:Ny, Nv)
    real :: Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv), Uxy(0:Nx, 0:Ny, Nv)
    real :: Sxx(0:Nx, 0:Ny), Syy(0:Nx, 0:Ny), Sxy(0:Nx, 0:Ny)
    real :: thick = 8d-3, Lx = 1.3, Ly = 3
    real :: x0 , xf , y0 , yf 
    integer :: i, j, i_pillar, j_pillar 
    integer, parameter :: Order = 4
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_pillar, y_pillar, R_pillar
    real :: dmin 
    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly
    x_pillar = Lx/2 
    y_pillar = 2*Ly/3 
    dmin = max( Lx/Nx, Ly/Ny) 
    R_pillar = max( 0.2, dmin ) 
    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1     
	 
       
!  call Non_Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
     
  call Grid_Initialization( "nonuniform", "x", Order, x )
  call Grid_Initialization( "nonuniform", "y", Order, y )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   
 
 Sxx =  -1d-6 * D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) ) * 6 / (thick)**2! + 1d-6 * ( Uyy(:,:,3))
 Syy =  -1d-6 * D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) ) * 6 / (thick)**2 !+ 1d-6 * ( Uxx(:,:,3))
 Sxy =  -1d-6 * D * (1 - poisson) * Uxy(:,:,1) * 6 / (thick)**2 !- 1d-6 * ( Uxy(:,:,3))

open(100, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Pillars_Non_Linear_Cartesian\Pillars_Non_Linear_Cartesian2.plt')
write (100,*) 'VARIABLES = "X", "Y", "Desplazamiento (m)", "sigma_xx", "sigma_yy", "sigma_xy" '
write (100,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (100,*) x(i), y(j), U(i,j,1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (100)

open(32, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Pillars_Non_Linear_Cartesian\Pillars_Non_Linear_Cartesian_phi.plt')
write (32,*) 'VARIABLES = "X", "Y", "phi_yy", "phi_xx", "phi_xy"'
write (32,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
 do  i=0 , Nx
  do j = 0, Ny
 write (32,*) x(i), y(j), 1d-6 *Uyy(i,j,3), 1d-6 *Uxx(i,j,3),- 1d-6 *Uxy(i,j,3)
  end do
 end do
close (32)
!****************** Uz PLOTTING **********************		

    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
    
    !****************** Uz PLOTTING **********************		
    

    
!    call scrmod('revers')
!    call metafl('xwin')
!    call disini() 
!    CALL TEXMOD('ON')
!    CALL TITLIN (" $  Sxx $ ", 2)
!	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)
	

contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy
  real :: phi, phixx, phiyy, phixy, F, Fxx, Fyy 
  real :: pillar_region , q
    
  pillar_region = (x - x_pillar)**2 +(y - y_pillar)**2 - R_pillar**2
  q = ( - rho_water * g * (y - y0))  
       
        v = u(1);   vxx = uxx(1);   vyy = uyy(1) ; vxy = uxy(1)
        w = u(2);   wxx = uxx(2);   wyy = uyy(2) 
      phi = u(3); phixx = uxx(3); phiyy = uyy(3) ; phixy = uxy(3)
        F = u(4);   Fxx = uxx(4);   Fyy = uyy(4) 
        
    if (pillar_region > 0) then   
        L(1) =  vxx +vyy - w   
        L(2) =  D*(wxx + wyy) - q - (thick) * ( phiyy*vxx + phixx*vyy - 2* phixy * vxy )
        L(3) =  phixx + phiyy - F
        L(4) =  Fxx + Fyy  + (E/2) *(thick**2) *  ( vyy*vxx + vxx*vyy - 2* vxy * vxy )
             
     else
        L(1) = vxx +vyy - w 
        L(2) = v 
        L(3) =  phixx + phiyy - F
        L(4) =  Fxx + Fyy  + (E/2) *(thick**2) *  ( vyy*vxx + vxx*vyy - 2* vxy * vxy )        
      
    endif  
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
            
        elseif (x==xf) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
            
        elseif (y==y0) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
        elseif (y==yf) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 


subroutine Test_Pillars_Polar_Navier_Analytical

    integer, parameter :: Nx = 40, Ny = 40, N_order = 15 
    real :: x(0:Nx), y(0:Ny), V(0:Nx), Vx(0:Nx), Vxx(0:Nx), Vxxx(0:Nx)
    real :: Mrr(0:Nx), Srr(0:Nx), Mtt(0:Nx), Stt(0:Nx), Qr(0:Nx), Srz(0:Nx)
    real :: thick = 8d-3, Array(5,0:Nx)
    real :: x0 = 0, xf = 0.5, y0, yf
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: h_pillar = 1.5, rho_water = 1000, p_atm = 101325 , g = 9.817, R_pillar = 10d-3
    real :: p_pillar , D , A, B, C, Kq, alpha
     yf = thick/2; y0 = -thick/2
     y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]
     x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
     call Grid_Initialization( "nonuniform", "x", N_Order, x )
    


    p_pillar = - rho_water * g * h_pillar

	D = E * (thick**3) / (12 * (1- poisson**2))
    alpha = (5 + poisson)/(3 + poisson)
    Kq = p_pillar * (xf**2)/(64*D)
    
    A = 2* Kq *alpha
    B = Kq*(2*alpha*log(xf)-1)
    do i = 0, Nx
      if ( x(i) > R_pillar ) then
     V(i) = (Kq / (xf**2)) * (x(i)- R_pillar)**4 - A * (x(i)- R_pillar)**2 * log((x(i)- R_pillar)) + B* (x(i)- R_pillar)**2
     
       else
       V(i) = 0
       end if
    end do
    

    call Derivative( "x", 1, V, Vx) 
    call Derivative( "x", 1, Vx, Vxx) 
    call Derivative( "x", 1, Vxx, Vxxx) 

  do i = 0, Nx
      if ( x(i) > R_pillar ) then
     Mrr(i) = -D*(Vxx(i) + poisson*Vx(i)/(x(i)-R_pillar))

     Mtt(i) = -D*(poisson*Vxx(i) + Vx(i)/(x(i)-R_pillar))
!~      Qr(i) = p_pillar * (x(i)-R_pillar) /2   - 4* A* D /(x(i)-R_pillar)   
     Qr(i) =  D*(Vxxx(i) + Vxx(i)/(x(i)-R_pillar) - Vx(i)/((x(i)-R_pillar)**2))
   
   
     else
     Mrr(i) =-D* (Vxx(i))
     Mtt(i) = -D*(poisson*Vxx(i))
     Qr(i) = D*(Vxxx(i) )    

      end if
  end do
Srr = 6 * Mrr /(thick**2)
Stt = 6 * Mtt /(thick**2)
Srz = Qr /(thick)

do i = 0,Nx
 Array(1,i) = x(i)
 Array(2,i) = V(i)
 Array(3,i) = Srr(i)
 Array(4,i) = Stt(i)
 Array(5,i) = Srz(i)
end do

!open(13, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\JAHR_Numerical - copy\Pillars\Graficas\Pillar_Polar_Navier_Analytical.plt')
!write (13,*) 'VARIABLES = "r", "Uz (m)", "\sigma_rr (Pa)", "\sigma_tt (Pa)", "\sigma_rz (Pa)",'
!write (13,*) 'ZONE DATAPACKING=POINT, I=', Nx
! do  i=0 , Nx
!  do j = 1, 5
! write (13,*) Array(j,i)
!  end do
! end do
!close (13) 

!------- Uz PLOTTING ----------------	
	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    CALL QPLOT (x, V*1d3, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')    
    CALL TITLIN (" $  Sigma RR max (MPa) $ ", 2)
    CALL QPLOT (x, Srr * 1d-6   , Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')    
    CALL TITLIN (" $  Sigma TT max (MPa) $ ", 2)
    CALL QPLOT (x, Stt * 1d-6   , Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')    
    CALL TITLIN (" $  Sigma RZ max (MPa) $ ", 2)
    CALL QPLOT (x, Srz * 1d-6   , Nx+1)
    


  

end subroutine  




subroutine Test_Pillars_Polar_system

    integer, parameter :: Nx = 55, Ny = 55, Nv = 2, N_order = 15 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    real :: Ur(0:Nx,0:Ny,Nv), Uz(0:Nx,0:Ny,Nv)
    real :: Srr(0:Nx, 0:Ny), Stt(0:Nx, 0:Ny), Szz(0:Nx, 0:Ny),Srz(0:Nx, 0:Ny)
    real :: Srr1D(0:Nx), Stt1D(0:Nx), Srz1D(0:Nx), U_z(0:Nx), U_zr(0:Nx), U_zrr(0:Nx)

    real :: thick = 8, D
    real :: x0 = 0, xf = 500, y0 , yf
    integer :: i, j, m
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, K
    real :: h_pillar = 1.5, rho_water = 1000, p_atm = 101325 , g = 9.817, R_pillar = 10
    real :: p_pillar  
    
    y0 = -0.5 * thick ; yf = 0.5 * thick
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]
    
    p_pillar = p_atm - rho_water * g * h_pillar

    mu = 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
	K = (mu) / (lambda + 2*mu)
    U = 1   
       
  !!  call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = N_order, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
      call Grid_Initialization( "nonuniform", "x", N_Order, x )
     call Grid_Initialization( "nonuniform", "y", N_Order, y )
   
    do m = 1, Nv
        call Derivative( ["x","y"], 1, 1, U(0:,0:,m), Ur(0:,0:,m)  )
        call Derivative( ["x","y"], 2, 1, U(0:,0:,m), Uz(0:,0:,m)  )
     end do
    
  do i = 0, Nx
    do j = 0, Ny 
    
     if (i>0) then
    
  Srr(i,j) = (lambda + 2*mu)* Ur(i,j,1) + lambda * (U(i,j,1)/x(i) + Uz(i,j,2))
  Stt(i,j) = (lambda + 2*mu)* U(i,j,1)/x(i) + lambda * (Ur(i,j,1) + Uz(i,j,2))
  Szz(i,j) = (lambda + 2*mu)* Uz(i,j,2) + lambda * (Ur(i,j,1) + U(i,j,1)/x(i))
  Srz(i,j) = 2*mu *(Ur(i,j,2) + Uz(i,j,1))
    
     else 
  Srr(i,j) = (lambda + 2*mu)* Ur(i,j,1) + lambda * ( Uz(i,j,2))
  Stt(i,j) =  lambda * (Ur(i,j,1) + Uz(i,j,2))
  Szz(i,j) = (lambda + 2*mu)* Uz(i,j,2) + lambda * (Ur(i,j,1) )
  Srz(i,j) = 2*mu *(Ur(i,j,2) + Uz(i,j,1))
  
     endif
     
    enddo
  enddo   
  D = E * (thick**3) / (12*(1-poisson**2))

  
!open(13, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\JAHR_Numerical - copy\Pillars\Graficas\Pillar_Polar_Elastic.plt')
!write (13,*) 'VARIABLES = "r", "z", "Ur (m)", "Uz (m)", "\sigma_rr (Pa)", "\sigma_tt (Pa)", "\sigma_zz (Pa)", "\sigma_rz (Pa)",'
!write (13,*) 'ZONE I=', Nx,', J=', Ny,',DATAPACKING=POINT'
! do  i=0 , Nx
!  do j = 0, Ny
! write (13,*) x(i), y(j), U(i,j,1), U(i,j,2), Srr(i,j),  Stt(i,j), Szz(i,j), Srz(i,j)
!  end do
! end do
!close (13)                   
!------- Uz PLOTTING ----------------	
  
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_Z (mm) $ ", 2)
    call qplcon( U(:,:,2), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_Z (mm) $ ", 2)
	CALL QPLCLR ( U(:,:,2), Nx+1, Ny+1)
	
	
!--------- Stress Srr---------	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rr} (Pa) $ ", 2)
	CALL qplcon ( Srr(:,:), Nx+1, Ny+1, 20)
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rr} (Pa) $ ", 2)
	CALL qplclr ( Srr(:,:), Nx+1, Ny+1)
	
!--------- Stress Stt---------	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{\theta \theta} (Pa) $ ", 2)
	CALL qplcon ( Stt(:,:), Nx+1, Ny+1, 20)
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{\theta \theta} (Pa) $ ", 2)
	CALL qplclr ( Stt(:,:), Nx+1, Ny+1)


!--------- Stress Szz---------	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{zz} (Pa) $ ", 2)
	CALL qplcon ( Szz(:,:), Nx+1, Ny+1, 20)
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{zz} (Pa) $ ", 2)
	CALL qplclr ( Szz(:,:), Nx+1, Ny+1)
	
!--------- Stress Srz---------	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rz}  $ ", 2)
	CALL qplcon ( Srz(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rz}  $ ", 2)
	CALL qplclr ( Srz(:,:), Nx+1, Ny+1)
	


!------- Ur PLOTTING ----------------	
	
	
	 call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_r (mm) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_r (mm) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	

contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
       real :: v, vx,vy, w, wx, wy
       real :: vxx, vyy , vxy
       real :: wxx, wyy, wxy
       
        v = u(1); vx = ux(1); vy = uy(1) 
        w = u(2); wx = ux(2); wy = uy(2)   
        
        vxx = uxx(1) ; vyy = uyy(1) ; vxy = uxy(1)      
        wxx = uxx(2) ; wyy = uyy(2) ; wxy = uxy(2)

 L(1) = (x**2)* vxx + (x**2)* K* vyy + (x**2)* (1-K)* wxy +  x* vx  -  v 
 L(2) = (x)* wyy + (x)*(K)* wxx + (x)* (1-K)* vxy + (1-K)* vy + K* wx 


        
end function

!********* Boundary conditions *********
function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = ux(1)
            BCs(2) = ux(2)
            
        elseif (x==xf) then  
          
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (y==y0) then
           if (x > R_pillar) then
          
            BCs(1) = uy(2) + (1-2*K) * ux(1) + (1-2*K) * u(1)/x + p_pillar /(lambda + 2*mu)
            BCs(2) = K * (uy(1) + ux(2)) 
           else 
            BCs(1) = u(1)
            BCs(2) = u(2) 
           endif
        elseif (y==yf) then
         if (x > R_pillar) then
  
            BCs(1) = x* uy(2) + x* (1-2*K)* ux(1) + (1-2*K)* u(1) + x* p_atm /(lambda + 2*mu)
            BCs(2) = K * (uy(1) + ux(2)) 
             else 
            BCs(1) = u(1)
            BCs(2) = u(2) 
           endif
            
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, yf=", x0, xf
            stop 
        endif

    end function

end subroutine  

subroutine Test_Pillars_Polar_Navier_system

    integer, parameter :: Nx = 40, Ny = 40, Nv = 2, Order = 15 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), V(0:Nx), Vx(0:Nx), Vxx(0:Nx) , Vxxx(0:Nx)
    real :: Mr(0:Nx), Mt(0:Nx), Qr(0:Nx), Srr(0:Nx), Stt(0:Nx), Srz(0:Nx)
    real :: thick = 8
    real :: x0 = 0, xf = 500, y0 , yf
    integer :: i, j, m
    real :: E = 72d9 , poisson = 0.22 
    real :: h_pillar = 1.5, rho_water = 1000, p_atm = 101325 , g = 9.817, R_pillar = 10
    real :: p_pillar , D 
    
    y0 = -0.5 * thick ; yf = 0.5 * thick
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]
    
    p_pillar = - rho_water * g * h_pillar

	D = E * (thick**3) / (12 * (1- poisson**2))

    U = 1      
  !  call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = N_order, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
  
     call Grid_Initialization( "nonuniform", "x", Order, x )
    V(0:Nx)   = U(0:Nx,0,1)
	Vxx(0:Nx) = U(0:Nx,0,2)
	
    call Derivative( "x", 1,   V,   Vx) 
    call Derivative( "x", 1, Vxx, Vxxx) 
!~ 
     
  
     
    do i = 0, Nx
      if (i >  0) then
     
     Mr(i) = -D *(Vxx(i) + poisson * Vx(i) / x(i))
     Mt(i) = -D *(poisson * Vxx(i) +  Vx(i) / x(i))
     Qr(i) =  D*(Vxxx(i) + Vxx(i)/x(i) - Vx(i)/(x(i)**2))
     
      else
      
     Mr(i) = -D *(Vxx(i))
     Mt(i) = -D *(poisson *  Vxx(i))  
     Qr(i) =  D*(Vxxx(i) )
  
      endif
    enddo
 !   Maximum stresses 
    Srr = 6 * Mr / (thick )**2
    Stt = 6 * Mt / (thick )**2
    Srz = Qr / (thick )
    
!~ open(13, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\JAHR_Numerical - copy\Pillars\Graficas\Pillar_Polar_Navier.plt')
!~ write (13,*) 'VARIABLES = "r", "z", "Uz (m)", "\sigma_rr (Pa)", "\sigma_tt (Pa)", "\sigma_rz (Pa)",'
!~ write (13,*) 'ZONE I=', Nx,', J=', Ny,',DATAPACKING=POINT'
!~  do  i=0 , Nx
!~   do j = 0, Ny
!~  write (13,*) x(i), y(j), U(i,j,1), Srr(i),  Stt(i), Srz(i)
!~   end do
!~  end do
!~ close (13) 
!------- Uz PLOTTING ----------------	
	
	
	 call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    CALL QPLOT (x, V, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rr} $ ", 2)
    CALL QPLOT (x, Srr, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{\theta \theta} $ ", 2)
    CALL QPLOT (x, Stt, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rz} $ ", 2)
    CALL QPLOT (x, Srz, Nx+1)
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
       real :: v, vx,vy, w, wx, wy
       real :: vxx, vyy , vxy
       real :: wxx, wyy, wxy
       
        v = u(1); vx = ux(1); vy = uy(1) 
        w = u(2); wx = ux(2); wy = uy(2)   
        vxx = uxx(1) ; vyy = uyy(1) ; vxy = uxy(1)      
        wxx = uxx(2) ; wyy = uyy(2) ; wxy = uxy(2)
        
   if (x > R_pillar) then
    L(1) = vxx - w
    L(2) = (x**3)* wxx + 2* x**2 * wx - x * vxx + vx - p_pillar * (x**3) /D !(0.65/1.65)*
   else
    L(1) = vxx-w
    L(2) = v
   endif
    

        
end function

!********* Boundary conditions *********
function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))
        
		if (x==x0) then
            BCs(1) = ux(1)
            BCs(2) = u(2)
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2) + poisson * ux(1) / xf
            
                    
        elseif (y==y0) then
      
            BCs(1) = uy(1)
            BCs(2) = uy(2) 
           
        elseif (y==yf) then
                        
            BCs(1) = uy(1)
            BCs(2) = uy(2)
            
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, yf=", x0, xf
            stop 
        endif

    end function

end subroutine  

subroutine Test_Pillars_Polar_Non_Linear_Navier_system

    integer, parameter :: Nx = 40, Ny = 40, Nv = 4, Order = 5 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), V(0:Nx), Vx(0:Nx), Vxx(0:Nx) , Vxxx(0:Nx)
    real :: Mr(0:Nx), Mt(0:Nx), Qr(0:Nx), Srr(0:Nx), Stt(0:Nx), Srz(0:Nx)
    real :: thick = 8
    real :: x0 = 0, xf = 500, y0 , yf
    integer :: i, j, m
    real :: E = 72d9 , poisson = 0.22 
    real :: h_pillar = 1.5, rho_water = 1000, p_atm = 101325 , g = 9.817, R_pillar = 20
    real :: p_pillar , D 
    
    y0 = -0.5 * thick ; yf = 0.5 * thick
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]
    
    p_pillar = - rho_water * g * h_pillar

	D = E * (thick**3) / (12 * (1- poisson**2))

    U = 0     
  !  call Non_Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = N_order, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
  
     call Grid_Initialization( "nonuniform", "x", Order, x )
    V(0:Nx)   = U(0:Nx,0,1)
	Vxx(0:Nx) = U(0:Nx,0,2)
	
    call Derivative( "x", 1,   V,   Vx) 
    call Derivative( "x", 1, Vxx, Vxxx) 
!~ 
     
  
     
    do i = 0, Nx
      if (i >  0) then
     
     Mr(i) = -D *(Vxx(i) + poisson * Vx(i) / x(i))
     Mt(i) = -D *(poisson * Vxx(i) +  Vx(i) / x(i))
     Qr(i) =  D*(Vxxx(i) + Vxx(i)/x(i) - Vx(i)/(x(i)**2))
     
      else
      
     Mr(i) = -D *(Vxx(i))
     Mt(i) = -D *(poisson *  Vxx(i))  
     Qr(i) =  D*(Vxxx(i) )
  
      endif
    enddo
 !   Maximum stresses 
    Srr = 6 * Mr / (thick )**2
    Stt = 6 * Mt / (thick )**2
    Srz = Qr / (thick )
    

!------- Uz PLOTTING ----------------	
	
	
	 call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    CALL QPLOT (x, V, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rr} $ ", 2)
    CALL QPLOT (x, Srr, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{\theta \theta} $ ", 2)
    CALL QPLOT (x, Stt, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rz} $ ", 2)
    CALL QPLOT (x, Srz, Nx+1)
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
       real :: v, vx,vy, w, wx, wy, phi, phix, phiy, F, Fx
       real :: vxx, vyy , vxy
       real :: wxx, wyy, wxy
       real :: phixx, phiyy, phixy, Fxx, Fyy
       
          v = u(1);   vx = ux(1);   vy = uy(1) 
          w = u(2);   wx = ux(2);   wy = uy(2)  
        phi = u(3); phix = ux(3); phiy = uy(3)   
          F = u(4) ;  Fxx = uxx(4)
        
          vxx = uxx(1) ;   vyy = uyy(1) ;   vxy = uxy(1)      
          wxx = uxx(2) ;   wyy = uyy(2) ;   wxy = uxy(2)
        phixx = uxx(3) ; phiyy = uyy(3) ; phixy = uxy(3)
          Fxx = uxx(4)
        
   if (x > R_pillar) then
    L(1) = vxx - w
    L(2) = (x**3)* wxx + 2* x**2 * wx - x * vxx + vx - p_pillar * (x**3) /D - (thick/D)*(  (x**2) * vxx * phix +  (x**2) *phixx * vx     ) 
    L(3) = phixx - F
    L(4) = (x**3)* Fxx + 2* x**2 * Fx - x * phixx + phix + (E/2)*(   (x**2) * vxx * vx +  (x**2) *vxx * vx            ) 
   else
    L(1) = vxx-w
    L(2) = v
    L(3) = phixx - F
    L(4) = (x**3)* Fxx + 2* x**2 * Fx - x * phixx + phix + (E/2)*(   (x**2) * vxx * vx +  (x**2) *vxx * vx            ) 

   endif   
        
end function

!********* Boundary conditions *********
function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))
        
		if (x==x0) then
            BCs(1) = ux(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2) + poisson * ux(1) / xf
            BCs(3) = u(3)
            BCs(4) = u(4)            
                    
        elseif (y==y0) then
      
            BCs(1) = uy(1)
            BCs(2) = uy(2) 
            BCs(3) = uy(3)
            BCs(4) = uy(4)        
        elseif (y==yf) then
                        
            BCs(1) = uy(1)
            BCs(2) = uy(2) 
            BCs(3) = uy(3)
            BCs(4) = uy(4)
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, yf=", x0, xf
            stop 
        endif

    end function

end subroutine  








subroutine Linear_Cylindrical_Plate_Validation

    integer, parameter :: Nx = 40, Ny = 40, Nv = 2, Order = 15 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), V(0:Nx), Vx(0:Nx), Vxx(0:Nx) , Vxxx(0:Nx)
    real :: Mr(0:Nx), Mt(0:Nx), Qr(0:Nx), Srr(0:Nx), Stt(0:Nx), Srz(0:Nx)
    real :: thick = 8
    real :: x0 = 0, xf = 500, y0 , yf
    integer :: i, j, m
    real :: E = 72d9 , poisson = 0.22 
    real :: h_pillar = 1.5, rho_water = 1000, p_atm = 101325 , g = 9.817, R_pillar = 10
    real :: p_pillar , D 
    
    y0 = -0.5 * thick ; yf = 0.5 * thick
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]
    
    p_pillar = - rho_water * g * h_pillar

	D = E * (thick**3) / (12 * (1- poisson**2))

    U = 1      
  !  call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = N_order, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
  
   call Grid_Initialization( "nonuniform", "x", Order, x )
    V(0:Nx)   = U(0:Nx,0,1)
	Vxx(0:Nx) = U(0:Nx,0,2)
	
    call Derivative( "x", 1,   V,   Vx) 
    call Derivative( "x", 1, Vxx, Vxxx) 
!~ 
     
  
     
    do i = 0, Nx
      if (i >  0) then
     
     Mr(i) = -D *(Vxx(i) + poisson * Vx(i) / x(i))
     Mt(i) = -D *(poisson * Vxx(i) +  Vx(i) / x(i))
     Qr(i) =  D*(Vxxx(i) + Vxx(i)/x(i) - Vx(i)/(x(i)**2))
     
      else
      
     Mr(i) = -D *(Vxx(i))
     Mt(i) = -D *(poisson *  Vxx(i))  
     Qr(i) =  D*(Vxxx(i) )
  
      endif
    enddo
 !   Maximum stresses 
    Srr = 6 * Mr / (thick )**2
    Stt = 6 * Mt / (thick )**2
    Srz = Qr / (thick )
    
!~ open(13, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\JAHR_Numerical - copy\Pillars\Graficas\Pillar_Polar_Navier.plt')
!~ write (13,*) 'VARIABLES = "r", "z", "Uz (m)", "\sigma_rr (Pa)", "\sigma_tt (Pa)", "\sigma_rz (Pa)",'
!~ write (13,*) 'ZONE I=', Nx,', J=', Ny,',DATAPACKING=POINT'
!~  do  i=0 , Nx
!~   do j = 0, Ny
!~  write (13,*) x(i), y(j), U(i,j,1), Srr(i),  Stt(i), Srz(i)
!~   end do
!~  end do
!~ close (13) 
!------- Uz PLOTTING ----------------	
	
	
	 call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (mm) $ ", 2)
    CALL QPLOT (x, V, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rr} $ ", 2)
    CALL QPLOT (x, Srr, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{\theta \theta} $ ", 2)
    CALL QPLOT (x, Stt, Nx+1)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  \sigma_{rz} $ ", 2)
    CALL QPLOT (x, Srz, Nx+1)
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
       real :: v, vx,vy, w, wx, wy
       real :: vxx, vyy , vxy
       real :: wxx, wyy, wxy
       
        v = u(1); vx = ux(1); vy = uy(1) 
        w = u(2); wx = ux(2); wy = uy(2)   
        vxx = uxx(1) ; vyy = uyy(1) ; vxy = uxy(1)      
        wxx = uxx(2) ; wyy = uyy(2) ; wxy = uxy(2)
        
   if (x > R_pillar) then
    L(1) = vxx - w
    L(2) = (x**3)* wxx + 2* x**2 * wx - x * vxx + vx - p_pillar * (x**3) /D !(0.65/1.65)*
   else
    L(1) = vxx-w
    L(2) = v
   endif
    

        
end function

!********* Boundary conditions *********
function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))
        
		if (x==x0) then
            BCs(1) = ux(1)
            BCs(2) = u(2)
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2) + poisson * ux(1) / xf
            
                    
        elseif (y==y0) then
      
            BCs(1) = uy(1)
            BCs(2) = uy(2) 
           
        elseif (y==yf) then
                        
            BCs(1) = uy(1)
            BCs(2) = uy(2)
            
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, yf=", x0, xf
            stop 
        endif

    end function

end subroutine  








subroutine Non_Linear_Cylindrical_Plate_Validation

    integer, parameter :: Nx = 15, Ny = 15, Nv = 4
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    real :: Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv), Uxy(0:Nx, 0:Ny, Nv)
    real :: Sxx(0:Nx, 0:Ny), Syy(0:Nx, 0:Ny), Sxy(0:Nx, 0:Ny)
    real :: thick = 16d-3, Lx = 1.3, Ly = 3
    real :: x0 , xf , y0 , yf 
    integer :: i, j, i_pillar, j_pillar 
    integer, parameter :: Order = 5
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_pillar, y_pillar, R_pillar
    real :: dmin 
    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly
    x_pillar = Lx/2 
    y_pillar = 2*Ly/3 
    dmin = max( Lx/Nx, Ly/Ny) 
    R_pillar = max( 0.2, dmin ) 
    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       
!  call Non_Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = Order, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
     
   call Grid_Initialization( "nonuniform", "x", Order, x )
   call Grid_Initialization( "nonuniform", "y", Order, y )
 
  

 
 !Sxx = - 1d-6 * D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) ) * 6 / (thick)**2 + 1d-6 * ( Uyy(:,:,3))
 !Syy = - 1d-6 * D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) ) * 6 / (thick)**2 + 1d-6 * ( Uxx(:,:,3))
 !Sxy = - 1d-6 * D * (1 - poisson) * Uxy(:,:,1) * 6 / (thick)**2 - 1d-6 * ( Uxy(:,:,3))

!open(13, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\JAHR_Numerical - copy\Pillars\Graficas\Pillar_Non_Linear_Navier.plt')
!write (13,*) 'VARIABLES = "X", "Y", "Desplazamiento (m)", "\sigma_xx (MPa)", "\sigma_yy (MPa)", "\sigma_xy (MPa)",'
!write (13,*) 'ZONE I=', Nx,', J=', Ny,',DATAPACKING=POINT'
! do  i=0 , Nx
!  do j = 0, Ny
! write (13,*) x(i), y(j), U(i,j,1), Sxx(i,j),  Syy(i,j), Sxy(i,j)
!  end do
! end do
!close (13)

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
    
    !****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,4), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,4), Nx+1, Ny+1)
	

contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy
  real :: phi, phixx, phiyy, phixy, F, Fxx, Fyy 
  real :: pillar_region , q
    
  pillar_region = (x - x_pillar)**2 +(y - y_pillar)**2 - R_pillar**2
  q = ( - rho_water * g * (y - y0))  
       
        v = u(1);   vxx = uxx(1);   vyy = uyy(1) ; vxy = uxy(1)
        w = u(2);   wxx = uxx(2);   wyy = uyy(2) 
      phi = u(3); phixx = uxx(3); phiyy = uyy(3) ; phixy = uxy(3)
        F = u(4);   Fxx = uxx(4);   Fyy = uyy(4) 
    if (pillar_region > 0) then   
        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D + (thick/D) * ( phiyy*vxx + phixx*vyy - 2* phixy * vxy )
        L(3) =  phixx + phiyy - F
        L(4) =  Fxx + Fyy  + (E/2) * ( vyy*vxx + vxx*vyy - 2* vxy * vxy )
        
        
     else
        L(1) = vxx +vyy - w 
        L(2) = v 
        L(3) =  phixx + phiyy - F
        L(4) =  Fxx + Fyy  + (E/2) * ( vyy*vxx + vxx*vyy - 2* vxy * vxy )        
      
    endif  
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
            
        elseif (x==xf) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
            
        elseif (y==y0) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
        elseif (y==yf) then
            
            BCs(1) = u(1)
            BCs(2) = u(2)
            BCs(3) = u(3)
            BCs(4) = u(4)
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 


subroutine Test_1_Stripe_Navier_system

    integer, parameter :: Nx = 20, Ny = 20, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 10d-3, Lx = 1.3, Ly = 3
    real :: x0 , xf , y0 , yf, h0 = 0.5
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_stripe, w_stripe

    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

     x_stripe = Lx/2
    w_stripe = Lx/Nx

    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

 !   call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 4, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
   
    call Grid_Initialization( "nonuniform", "x", 4, x )
    call Grid_Initialization( "nonuniform", "y", 4, y )
 
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      


 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\Stripes_Belen\Stripe1_05.plt')
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, vx
  real :: q, stripe_region
    
  stripe_region = (x-x_stripe)**2 -(w_stripe**2)

  q = ( - rho_water * g * (y - h0))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2)  ; vx = ux(1)
       
 !   if ( (x_stripe + 0.5*Lx/Nx > x >= x_stripe) .or. (x_stripe - 0.5*Lx/Nx < x <= x_stripe) ) then   
    if ( stripe_region < 0 ) then   
     
        L(1) =  vxx +vyy - w   

        L(2) = v 
       
     else

        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
    endif  
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 



subroutine Test_2_Stripe_Navier_system

    integer, parameter :: Nx = 20, Ny = 20, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 10d-3, Lx = 1.3, Ly = 3
    real :: x0 , xf , y0 , yf, h0 = 0.5
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_stripe1, x_stripe2, w_stripe

    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

     x_stripe1 = Lx/3
     x_stripe2 = 2*Lx/3

    w_stripe = Lx/Nx

    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

!    call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 4, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
   
  call Grid_Initialization( "nonuniform", "x", 4, x ) 
  call Grid_Initialization( "nonuniform", "y", 4, y ) 
!  call Check_grid( "x", x, 4, size( U(:,0,1) ) ) 
!  call Check_grid( "y", y, 4, size( U(0,:,1) ) )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      


 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\Stripes_Belen\Stripe2_05.plt') 
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, vx
  real :: q, stripe_region1, stripe_region2
    
  stripe_region1 = (x-x_stripe1)**2 -(w_stripe**2)
  stripe_region2 = (x-x_stripe2)**2 -(w_stripe**2)


  q = ( - rho_water * g * (y - h0))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2)  ; vx = ux(1)
       
 !   if ( (x_stripe + 0.5*Lx/Nx > x >= x_stripe) .or. (x_stripe - 0.5*Lx/Nx < x <= x_stripe) ) then   
    if ( stripe_region1 < 0 ) then   
     
        L(1) =  vxx +vyy - w   

        L(2) = v 

    elseif ( stripe_region2 < 0 ) then   
     
        L(1) =  vxx +vyy - w   

        L(2) = v 
       
     else

        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
    endif  
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 



subroutine Test_Plate_05_Navier_system

    integer, parameter :: Nx = 20, Ny = 20, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 10d-3, Lx = 0.65, Ly = 3, h0 = 0.5
    real :: x0 , xf , y0 , yf
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_stripe, w_stripe

    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

  !  call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 4, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
   
  call Grid_Initialization( "nonuniform", "x", 4, x ) 
  call Grid_Initialization( "nonuniform", "y", 4, y ) 
  
 ! call Check_grid( "x", x, 4, size( U(:,0,1) ) ) 
 ! call Check_grid( "y", y, 4, size( U(0,:,1) ) )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      


 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\Stripes_Belen\Stripe05_05.plt')
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, vx
  real :: q, stripe_region
    

        q = ( - rho_water * g * (y - h0))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2)  ; vx = ux(1)
       

        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
        

end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 

subroutine Test_Plate_033_Navier_system

    integer, parameter :: Nx = 22, Ny = 22, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 10d-3, Lx = 0.43333, Ly = 3
    real :: x0 , xf , y0 , yf, h0 =0.5
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_stripe, w_stripe

    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

 !   call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 4, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
   
  call Grid_Initialization( "nonuniform", "x", 4, x ) 
  call Grid_Initialization( "nonuniform", "y", 4, y )  
  
!  call Check_grid( "x", x, 4, size( U(:,0,1) ) ) 
!  call Check_grid( "y", y, 4, size( U(0,:,1) ) )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      


 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\Stripes_Belen\Stripe033_05.plt')
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, vx
  real :: q, stripe_region
    

        q = ( - rho_water * g * (y - h0))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2)  ; vx = ux(1)
       

        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
        

end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)   
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = u(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

end function

end subroutine 



subroutine Test_Plate_Clamped_05_Navier_system

    integer, parameter :: Nx = 36, Ny = 36, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 10d-3, Lx = 0.65, Ly = 3
    real :: x0 , xf , y0 , yf
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_stripe, w_stripe

    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

   ! call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 6, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
   
  call Grid_Initialization( "nonuniform", "x", 6, x ) 
  call Grid_Initialization( "nonuniform", "y", 6, y )
  
!  call Check_grid( "x", x, 6, size( U(:,0,1) ) ) 
 ! call Check_grid( "y", y, 6, size( U(0,:,1) ) )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      


 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Stripes_Navier_Cartesian\Plate0.5_Clamped_Navier_Cartesian.plt')
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, vx
  real :: q, stripe_region
    

        q = ( - rho_water * g * (y - Ly*0.5))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2)  ; vx = ux(1)
       
    if (xf <= x + Lx/Nx) then

        L(1) =  vxx + vyy - w   
        L(2) =  v
     else   
        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
     end if
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) = u(1)
            BCs(2) = ux(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 

subroutine Test_Plate_Clamped_033_Navier_system

    integer, parameter :: Nx = 36, Ny = 36, Nv = 2 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv)
    real :: Uy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv)
    real :: Uxy(0:Nx, 0:Ny, Nv)
    real :: Mxx(0:Nx,0:Ny), Myy(0:Nx,0:Ny), Mxy(0:Nx,0:Ny)
    real :: Sxx(0:Nx,0:Ny), Syy(0:Nx,0:Ny), Sxy(0:Nx,0:Ny)
    real :: thick = 7d-3, Lx = 1.3/3, Ly = 3
    real :: x0 , xf , y0 , yf
    integer :: i, j
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D
    real :: rho_water = 1000, p_atm = 101325 , g = 9.817
    real :: x_stripe, w_stripe

    
    
    x0 = 0 ; xf =  Lx
    y0 = 0 ; yf =  Ly

    
    
    !Malla equiespaciada
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    U = 1      
	 
       

  !  call Linear_Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 4, N_variables = Nv, Differential_operator = L, Boundary_conditions = BCs, Solution = U )
     
  call Grid_Initialization( "nonuniform", "x", 4, x ) 
  call Grid_Initialization( "nonuniform", "y", 4, y )  
  
 ! call Check_grid( "x", x, 4, size( U(:,0,1) ) ) 
 ! call Check_grid( "y", y, 4, size( U(0,:,1) ) )
  
     do i=1, Nv 
        call Derivative( ["x","y"], 1, 1, U(0:,0:,i), Ux(0:,0:,i)  )
        call Derivative( ["x","y"], 1, 2, U(0:,0:,i), Uxx(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, U(0:,0:,i), Uy(0:,0:,i)  )
        call Derivative( ["x","y"], 2, 2, U(0:,0:,i), Uyy(0:,0:,i) )

        call Derivative( ["x","y"], 2, 1, Ux(0:,0:,i), Uxy(0:,0:,i))
     end do   

 Mxx = - D * (Uxx(:,:,1) + poisson * Uyy(:,:,1) )
 Myy = - D * (Uyy(:,:,1) + poisson * Uxx(:,:,1) )
 Mxy = - D * (1 - poisson) * Uxy(:,:,1)


do i = 0,Nx
    do j = 0,Ny
 Sxx(i,j) = 1d-6 *( Mxx(i,j) )* 6 / (thick)**2
 Syy(i,j) = 1d-6 *( Myy(i,j) )* 6 / (thick)**2
 Sxy(i,j) = 1d-6 *( Mxy(i,j) )* 6 / (thick)**2
    end do
end do

!****************** Uz PLOTTING **********************		
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
    call qplcon( U(:,:,1), Nx+1, Ny+1, 20)
    
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  Desplazamiento U_z (m) $ ", 2)
	CALL QPLCLR ( U(:,:,1), Nx+1, Ny+1)
	
!******************************************************	
!**************   Mxx  ********************************
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xx} $ ", 2)
	CALL qplcon ( Sxx(:,:), Nx+1, Ny+1, 20)
	
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xx} $ ", 2)
	CALL QPLCLR ( Sxx(:,:), Nx+1, Ny+1)

!******************************************************	
!**************   Myy  ********************************	
	
	call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{yy} $ ", 2)
	CALL qplcon ( Syy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{yy} $ ", 2)
	CALL QPLCLR ( Syy(:,:), Nx+1, Ny+1)
!******************************************************	
!**************   Mxy  ********************************
call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  S_{xy} $ ", 2)
	CALL qplcon ( Sxy(:,:), Nx+1, Ny+1, 20)

    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $ S_{xy} $ ", 2)
	CALL QPLCLR ( Sxy(:,:), Nx+1, Ny+1)

      


 open(101, file = 'C:\Users\Javi\Dropbox\Javier_Escoto\jahr\TFG-Francisco Javier Escoto López\figuras\Graficas\Stripes_Navier_Cartesian\Plate0.33_Clamped_Navier_Cartesian.plt')
write (101,*) 'VARIABLES = "X", "Y","U", "Sxx", "Syy", "Sxy" '
write (101,*) 'ZONE I=', Nx+1,', J=', Ny+1,',DATAPACKING=POINT'
  do i = 0, Nx 
  do  j =0 , Ny
 write (101,*) x(i), y(j), U(i,j, 1), Sxx(i,j), Syy(i,j), Sxy(i,j)
  end do
 end do
close (101) 
    
contains

!********* Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, vx
  real :: q, stripe_region
    

        q = ( - rho_water * g * (y - Ly*0.5))  
       
        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2)  ; vx = ux(1)
       

        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q / D
        

end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function BCs(x, y, u, ux, uy)

    
        real, intent(in) :: x, y, u(:), ux(:), uy(:)
        real :: BCs(size(u))

        if (x==x0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
            
        elseif (x==xf) then
            BCs(1) =  u(1)
            BCs(2) = ux(2)

        elseif (y==y0) then
            BCs(1) = u(1)
            BCs(2) = u(2)
      
        elseif (y==yf) then
            BCs(1) = u(1)
            BCs(2) = u(2)
     
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " x0, xf=", x0, xf
            write(*,*) " y0, yf=", y0, yf

            stop 
        endif

    end function

end subroutine 






!*******************************************************************************************************************************************
! Linear_1D
!*******************************************************************************************************************************************  
!subroutine Legendre_equation_1DP
!
!    integer, parameter :: N = 30 
!    real :: x(0:N), U(0:N) 
!    real :: x0 = -1 , xf = 1
!    integer :: i
!    real :: pi = 4 * atan(1.0)  
!
!    x = [ (x0 + (xf-x0)*i/N, i=0, N) ]
!
!    call Linear_Boundary_Value_Problem1DP( x_nodes = x, Order = 4,    & 
!                                        Differential_operator = L, & 
!                                        Boundary_conditions = BCs, & 
!                                        Solution = U )
!    call scrmod("reverse")
!    call qplot(x, U, N+1)
!
!contains 
!
!
!
!
!
!
!!line 50****** Differential operator *********
!    real function L(x, y, yx, yxx) 
!    
!        real, intent(in) :: x, y, yx, yxx                 
!        real, parameter :: n = 3.
!           
!      ! Legendre differential equation
!        L = (1. - x**2) * yxx - 2 * x * yxx + n * (n + 1.) * y
!           
!    end function 
!    
!!********* Boundary conditions *********
!    real function BCs(x, y, yx) 
!    
!        real, intent(in) :: x, y, yx            
!
!        if (x==x0) then
!                           BCs = y + 1
!        elseif (x==xf) then
!                           BCs = y - 1
!        else 
!            write(*,*) " Error BCs x=", x  
!            write(*,*) " a, b=", x0, xf
!            stop  
!        endif            
!                 
!    end function  
!
!end subroutine 



end module
