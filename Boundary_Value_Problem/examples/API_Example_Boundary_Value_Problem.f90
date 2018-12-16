module API_Example_Boundary_Value_Problem

    use Boundary_value_problems 
    use Finite_differences
    
    implicit none
    
    contains  
  
    
subroutine BVP_examples 

   call Legendre_1D
   call Linear_BVP2D
   
   call Linear_Plate_2D
   call Non_Linear_Plate_2D

end subroutine


!*******************************************************************************************************************************************
!  line 24 
!*******************************************************************************************************************************************  
subroutine Legendre_1D

    integer, parameter :: N = 30, q = 6  
    real :: x(0:N), U(0:N) 
    real :: x0 = -1 , xf = 1
    integer :: i
    real :: pi = 4 * atan(1.0)  

    x(0) = x0; x(N) = xf  
    call Grid_Initialization( grid_spacing = "nonuniform", &
                                 direction = "x",   q = q, nodes = x )
   
    call Boundary_Value_Problem( x_nodes = x, Order = q,               & 
                                 Differential_operator = Legendre,     & 
                                 Boundary_conditions   = Legendre_BCs, & 
                                 Solution = U )
    call scrmod("reverse")
    call qplot(x, U, N+1)

contains 





!line 50****** Differential operator *********
real function Legendre(x, y, yx, yxx) result(L)

    real, intent(in) :: x, y, yx, yxx                 
    real, parameter :: n = 3.
       
  ! Legendre differential equation
    L = (1. - x**2) * yxx - 2 * x * yx + n * (n + 1.) * y
       
end function 
    
!********* Boundary conditions *********
    real function Legendre_BCs(x, y, yx) result(BCs)
    
        real, intent(in) :: x, y, yx            

        if (x==x0) then
                           BCs = y + 1
        elseif (x==xf) then
                           BCs = y - 1
        else 
            write(*,*) " Error BCs x=", x  
            write(*,*) " a, b=", x0, xf
            stop  
        endif            
                 
    end function  

end subroutine 

!*******************************************************************************************************************************************
! Linear_2D
!*******************************************************************************************************************************************
subroutine Linear_BVP2D

    integer, parameter :: Nx = 30, Ny = 30, q = 11
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny)
    integer :: i, j
    real :: a=0, b=1
    real :: PI = 4 * atan(1.0)

    x(0) = a; x(Nx) = b; y(0) = a; y(Nx) = b; 
    
    call Grid_Initialization( "nonuniform", "x", q, x )
    call Grid_Initialization( "nonuniform", "y", q, y )
     
    U = 1      

    call Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = 11,  & 
                                 Differential_operator = L, Boundary_conditions = BCs, Solution = U ) 
                          
   
    call scrmod('revers')
    call metafl('xwin')
    call disini() 
    CALL TEXMOD('ON')
    CALL TITLIN (" $  u_{xx} + u_{yy}  = 50  \quad sin(2 \pi x) \quad sin( 2 \pi y )  $ ", 2)
    CALL TITLIN (" $ u(0,y) = 0, \quad u(1,y) = 0, \quad u(x,0) = 0, \quad u(x,1) = 0 $  ", 4)
    call qplcon( U, Nx+1, Ny+1, 20)

contains

!********* Function *********
    real function L(x, y, u, ux, uy, uxx, uyy, uxy)
    
        real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
        
        
        L =  uxx + uyy - 50 * sin(2*PI*x) * sin(2*PI*y) 

    end function

!********* Boundary conditions *********
    real function BCs(x, y, u, ux, uy)
    
        real, intent(in) :: x, y, u, ux, uy

        if (x==a) then
                          BCs = u
        elseif (x==b) then
                          BCs = u 
        elseif (y==a) then
                          BCs = u
        elseif (y==b) then
                          BCs = u 
        else
            write(*,*) " Error BCs x=", x
            write(*,*) " a, b=", a, b
            stop 
        endif

    end function

end subroutine 






!** line 1711
subroutine Linear_Plate_2D

    integer, parameter :: Nx = 20, Ny = 20, Nv = 2, q= 4  
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    real :: x0 = -1 , xf = 1 , y0 = -1  , yf = 1,  PI = 4* atan (1d0)
    integer :: i, j
  
    x(0) = x0 ; x(Nx) =  xf; y(0) = y0; y(Ny) = yf 
    call Grid_Initialization( "nonuniform", "x", q, x )
    call Grid_Initialization( "nonuniform", "y", q, y )
        
    U = 1          

    call Boundary_Value_Problem( x_nodes = x, y_nodes = y, Order = q,  & 
                                 N_variables = Nv,                     &
                                 Differential_operator = Linear_Plate, & 
                                 Boundary_conditions = L_Plate_BCs,    &
                                 Solution = U )
  
contains














!line 190***** Function *********
function Linear_Plate(x, y, u, ux, uy, uxx, uyy, uxy) result(L)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, q


        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2) 
       
        q =  4*pi**4 * sin(pi*x) * sin(pi*y) 

        L(1) =  vxx + vyy - w   
        L(2) =  wxx + wyy - q 
       

      
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function L_Plate_BCs(x, y, u, ux, uy) result(BCs)
    
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
            stop 
        endif

end function

end subroutine 

! line 242
subroutine Non_Linear_Plate_2D

    integer, parameter :: Nx = 20, Ny = 20, Nv = 4
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    real :: thick = 16d-3 
    real :: x0 = -1 , xf = 1 , y0 = -1 , yf = 1 
    integer :: i, j
    integer :: q = 6 ! Order 
    real :: E = 72d9 , poisson = 0.22 
    real :: mu, lambda, D, pi

    pi = acos(-1d0)
    
    
     
    
     x(0) = x0 ; x(Nx) =  xf; y(0) = y0; y(Ny) = yf 
     call Grid_Initialization( "nonuniform", "x", q, x )
     call Grid_Initialization( "nonuniform", "y", q, y )

    mu= 0.5 * E / (1+poisson) 
    lambda = E * poisson / ((1 - 2*poisson)*(1 + poisson))
    D = E * (thick**3) / (12 * (1- poisson**2))
    
    U=0
    
    linear2D = .false.  !TO eliminate ; do it automatically 
    
    call Boundary_Value_Problem( x_nodes = x, y_nodes = y,             & 
                                 Order = q,                            & 
                                 N_variables = Nv,                     & 
                                 Differential_operator = NL_Plate,     & 
                                 Boundary_conditions   = NL_Plate_BCs, & 
                                 Solution = U )
   
    call scrmod("reverse")
    call qplcon( U, Nx+1, Ny+1, 20) 
    
contains


!*line 284****** Function *********
function NL_Plate(x, y, u, ux, uy, uxx, uyy, uxy) result(L)

  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real ::   v,   vxx,   vyy,   vxy, w, wxx, wyy
  real :: phi, phixx, phiyy, phixy, F, Fxx, Fyy, q
  
  q = 1d-3 *D*((4*pi**4)) * sin(pi*x) * sin(pi*y) 
      
  v = u(1);   vxx = uxx(1);   vyy = uyy(1) ;   vxy = uxy(1)
  w = u(2);   wxx = uxx(2);   wyy = uyy(2) 
  phi = u(3); phixx = uxx(3); phiyy = uyy(3) ; phixy = uxy(3)
  F = u(4);   Fxx = uxx(4);   Fyy = uyy(4) 
  
  L(1) =  vxx + vyy - w   
  L(2) =  wxx + wyy - q/D + (thick/D)*(phiyy*vxx + phixx*vyy - 2*phixy*vxy)                                
  L(3) =  phixx + phiyy - F
  L(4) =    Fxx + Fyy  + (E/2) * ( vyy*vxx + vxx*vyy - 2* vxy * vxy )
 
end function

        
!**line 309**** Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function NL_Plate_BCs(x, y, u, ux, uy) result(BCs)

    
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




end module 












