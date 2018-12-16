module API_Example_Validation_BVP
    
    use Boundary_value_problems 
   ! use Finite_differences
    use BVP_Validation
        
implicit none

    contains
    
subroutine BVP_validation_examples

   call SHM_1D

   !call Linear_BVP2D
   
   !call Linear_BVP2D_System
   

end subroutine




subroutine SHM_1D

    real :: x(0:1)
    real, allocatable :: log_e(:), log_n(:)
    real :: x0 = 0 , xf = 2
    integer :: q, i
    real :: pi = 4 * atan(1.0)  
    

    x(0) = x0; x(1) = xf  

    do q=2, 8, 2
        call BVP_1D_Validation( Differential_operator = L, &
                                Boundary_conditions = BCs, &
                                Order = q, Spatial_Domain = x, &
                                Analytical_Solution = Analytical_Solution,&
                                log_e = log_e, log_N = log_N ) 
    
    
        call Save_validation_1D_linear(log_N, log_e, size(log_N), q)
    end do
contains 

!*****     LATEX - 50    *****
!***** Differential operator *********
    real function L(x, y, yx, yxx) 
    
        real, intent(in) :: x, y, yx, yxx 
     
      ! SHM differential equation
        L = yxx + (pi/4)**2 * y
           
    end function 
    
!*****     LATEX - 61     *****
!********* Boundary conditions *********
    real function BCs(x, y, yx) 
    
        real, intent(in) :: x, y, yx            

        if (x==x0) then
                           BCs = y 
        elseif (x==xf) then
                           BCs = y - 1
        else 
            write(*,*) " Error BCs x=", x  
            write(*,*) " a, b=", x0, xf
            stop  
        endif            
                 
    end function 
    
!*****     LATEX - 79     *****
    function Analytical_Solution(x) result(F)
    
        real, intent(in) :: x(:)
        real :: F( size(x) )
        
        integer :: i
        
        do i=1, size(x)
            F(i) =sin( pi/4*x(i) )
        end do
    
    end function

end subroutine


subroutine Linear_BVP2D

    integer :: q
    real :: x(0:1), y(0:1)
    real, allocatable :: log_e(:), log_n(:)
    
    integer :: i, j
    real :: a=0, b=1
    real :: PI = 4 * atan(1.0)

    x(0) = a; x(1) = b; y(0) = a; y(1) = b; 
    
    do q=2, 8, 2
        write(*,*) 'Order = ', q
        call BVP_2D_Validation( Differential_operator = L, Boundary_conditions = BCs, Order = q,                            &
                                Spatial_Domain_x = x, Spatial_Domain_y = y, Analytical_Solution = Analytical_Solution_2D,   &
                                log_e = log_e, log_N = log_N ) 
    
        write(*,*) 'log N = ', log_N
        write(*,*) 'log e = ', log_e
        call Save_validation_2D(log_N, log_e, size(log_N), q)
    end do
     


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
    
    function Analytical_Solution_2D(x, y) result(F)
    
        real, intent(in) :: x(:), y(:)
        real :: F( size(X), size(y) )
        
        integer :: i, j
        
        do i=1, size(x)
            do j=1, size(y)
                F(i,j) =- 50 / (2*(2*pi)**2) * sin( 2*pi*x(i) ) * sin( 2*pi*y(j) )
            end do
        end do
    
    end function

end subroutine 










!*****     LATEX - 180     *****
subroutine Linear_BVP2D_System

    integer, parameter ::  Nv = 2
    integer :: q
    real :: x(0:1), y(0:1)
    real, allocatable :: log_e(:), log_n(:)
    
    real :: x0 = -1 , xf = 1 , y0 = -1  , yf = 1 
    integer :: i, j
    real :: PI = 4* atan (1d0)
  
    x(0) = x0 ; x(1) =  xf; y(0) = y0; y(1) = yf 

     
    do q=2, 8, 2
        
        write(*,*) 'Order = ', q
        call BVP_2D_System_Validation( Differential_operator = L, Boundary_conditions = BCs, Order = q, Nvariables = Nv, Spatial_Domain_x = x, Spatial_Domain_y = y, Analytical_Solution = Analytical_Solution_2D_system, log_e = log_e, log_N = log_N ) 
    
        write(*,*) 'N = ', log_N
        write(*,*) 'e = ', log_e
        call Save_validation_2D_system(log_N, log_e, size(log_N), q)
        
    end do

  
    
    
contains
!*****     LATEX - 210     *****
!***** Function *********
function L(x, y, u, ux, uy, uxx, uyy, uxy)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u)) 
  real :: v, vxx, vyy, vxy, w, wxx, wyy, q


        v = u(1); vxx = uxx(1); vyy = uyy(1) ; vxy = uxy(1)
        w = u(2); wxx = uxx(2); wyy = uyy(2) 
       
        q =  4*pi**4 * sin(pi*x) * sin(pi*y) 

        L(1) =  vxx +vyy - w   
        L(2) =  wxx + wyy - q 
    
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
            stop 
        endif

end function

function Analytical_Solution_2D_System(x, y) result(F)
    
        real, intent(in) :: x(:), y(:)
        real :: F( size(X), size(y) )
        
        integer :: i, j
        
        do i=1, size(x)
            do j=1, size(y)
                F(i,j) = sin( pi*x(i) ) * sin( pi*y(j) )
            end do
        end do
    
end function
end subroutine 


end module
